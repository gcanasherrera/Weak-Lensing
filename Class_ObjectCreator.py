# Name: Class_ObjectCreator.py
#
# Weak-Lensing Validation Program III
#
# Type: Python class
#
# Content: 2 Classes, 1 constructor, 16 methods
#
# Description: General Class destinated to produce a simulation over a .fits picture in order to determinate whether Source Extractor detects properly the celestial objects.
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "2.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"


import matplotlib.pyplot as plt #Plot Libraries
import numpy as np #Maths arrays and more, matlab-type vectors/arrays
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
from astropy.modeling import models, fitting #Package for fitting functions with a astronomical character
import warnings #Advices
from astropy.io import fits #Open and Reading FITS Files usign astropy
import random #pseudo-random generatorCcCl
import seaborn as sns #Improvements for statistical-plots
from pymodelfit import FunctionModel1DAuto #Create own model of fitting
from scipy import spatial #KDTREE altorithm
import matplotlib.pylab as P #histograms


"""
     Class that defines the relation between F(m) according to the astrophysical definition
"""
class MagnitudeExponential(FunctionModel1DAuto):
    def f(self,m,a=1,b=5):
        return a*10**(-m/2.5 +b)


"""
General Class that contents several simulation methods destinated to introduce new celestial objects in a .fits picture with an ellipticial shape. There are different posibilities of simulations.
    
"""
class ObjectCreator(object):
    
    
    
    """
        Constructor that defines attribute of future class-objects

    """
    
    def __init__(self, fcat):  #define attributes of future class-objects
        
        self.fcat = fcat
        self.simulation_picture = ''
    
        self.mean_mag = np.mean(fcat['mag_iso'])
        self.error_mean_mag =np.std(fcat['mag_iso'])/math.sqrt(len(fcat['mag_iso']))
        self.mean_intensity = 0.0
        self.mean_ellip = np.mean(fcat['ellipticity'])
        self.mean_b = np.mean(fcat['B'])
        self.mean_ec = np.sqrt(1-(1-self.mean_ellip)*(1-self.mean_ellip))
        self.mean_a = self.mean_b/math.sqrt(1-self.mean_ec*self.mean_ec)
        
        self.error_mean_b = np.std(fcat['B'])/math.sqrt(len(fcat['B']))
        self.error_mean_ec = self.mean_ec*0.5*2*1/((1-self.mean_ec*self.mean_ec)*(1-self.mean_ec*self.mean_ec)*(1-self.mean_ec*self.mean_ec))*math.sqrt((1-self.mean_ec*self.mean_ec)*(1-self.mean_ec*self.mean_ec))
                            
        self.error_mean_a = self.mean_a*math.sqrt((self.error_mean_b/self.mean_b)*(self.error_mean_b/self.mean_b)+(self.mean_ec/self.error_mean_ec)*(self.mean_ec/self.error_mean_ec))
        
        self.parameter_a = 0.0
        self.parameter_b = 0.0


        self.out_mag_all = []
        self.out_mag = []
        self.out_flux = []
        self.out_flux_max = []
        self.out_mag_after_transf = []
        
        
        self.posible_obj_distances = []
        self.posible_obj_index = []
        self.cont_lost_obj_per = []


        self.random_mag = []
        #self.random_intensity = []
        #self.random_ellip = []
        #self.random_b = []
        #self.random_a = []
    
    
        self.x_position_int = []
        self.y_position_int = []
        self.x_position_simulation = []
        self.y_position_simulation = []
    
        self.x_data_image = 0
        self.y_data_image = 0
        self.picture_data = []
        self.matrix_data = []
        
        
        self.total = 0
        self.cont_1 = 0
        self.cont_0 = 0
        self.percentage_1 = 0
        self.percentage_0 = 0
        self.packing = 0.0
        self.number_to_packing = 0
        self.lost_objects = []
    
    
    """
        Method that fits the magnitude as a function of the flux in logaritmic scale through a third-order polynomical model
        
    """

    def plot_fitting_pol(self, mag_1, mag_2, val_1, val_2, lin):
        plt.figure()
        sns.set_style("white")
        sns.set_style("ticks")
        plt.xlabel(('log({})').format(val_1))
        plt.ylabel(('log({})').format(val_2))
        plt.title(('{} Vs. {} in Logaritmic Scale').format(val_1, val_2))
        xp = np.linspace(1, lin, 100)
        fit = np.polyfit(np.log(mag_1), np.log(mag_2), 3)
        p3 = np.poly1d(fit)
        print 'The fit of {} Vs. {} is: {}\n'.format(val_1, val_2, p3)
        plt.plot(np.log(mag_1), np.log(mag_2), 'r--', label='Source Extractor')
        #plt.errorbar(np.log(fcat['flux_iso']), np.log(fcat['mag_iso']), xerr = xerr_norm, yerr = yerr_norm, fmt='r.')
        plt.plot(xp, p3(xp), 'k-', label='fitting')
        plt.legend()
        plt.show(block=False)
        return p3
    
    """
        Method that fits the magnitude as a function of the flux in normal scale through an exponential model defined at http://astropy.readthedocs.org/en/latest/api/astropy.modeling.powerlaws.ExponentialCutoffPowerLaw1D.html#astropy.modeling.powerlaws.ExponentialCutoffPowerLaw1D.param_names
        
    """

    def plot_fitting_exp(self, x, y):
        sns.set(style="white", palette="muted", color_codes=True)
        t_init = models.ExponentialCutoffPowerLaw1D(amplitude=10., x_0=10., alpha=100., x_cutoff=10)
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, x, y)
        plt.figure()
        
        plt.loglog(x, y, 'ko', label='Source Extractor')
    
        plt.errorbar(x, y, self.fcat['fluxerr_iso'], self.fcat['magger_iso'], fmt='ko', label='Source Extractor')
        #plt.plot(x, y, 'ko', label='Source Extractor')
        xp= np.linspace(8.2, 25, len(x))
        plt.loglog(xp, t(xp), 'r--', label='Exponential Fitting')
        plt.xlabel('mag(iso)')
        plt.ylabel('flux(iso)')
        plt.legend()
        plt.title('Mag_iso Vs. Flux_iso')
        plt.show()
        print 'The fit of Flux_iso Vs. mag_iso is: {}\n'.format(t)
        #print t.alpha
        self.mean_intensity = t(self.mean_mag)
    #print 'The mean intensity is : {}\nThe mean isophotal magnitude is : {} '.format(self.mean_intensity, self.mean_mag)
    
    
    
    """
        Method that fits the magnitude as a function of the flux in normal scale through an exponential model defined myself using
        
    """
    
    def plot_fitting_exp_owndefined(self, x, y):
        sns.set(style="white", palette="muted", color_codes=True)
        model = MagnitudeExponential()
        af,bf = model.fitData(x,y)
        model.plot()
        #plt.errorbar(x, y, self.fcat['magger_iso'], self.fcat['fluxerr_iso'], fmt='r.')

        plt.xlabel('$m(iso)$')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylabel('$F(iso)$')
        
        plt.show()
        self.mean_intensity = model.f(m = self.mean_mag, a = af, b = bf)
        print 'The mean intensity is : {}\nThe mean isophotal magnitude is : {} '.format(self.mean_intensity, self.mean_mag)
        print 'The error in mean intensity is : {}\nThe error in mean isophotal magnitude is : {} '.format(self.error_mean_mag, self.error_mean_mag)
        
        
        return model, af, bf
    
    
    
    """
        Method that plots the histograms corresponding to the magnitudes: flux_iso, mag_iso, ellipticity, A (major axis of the ellipse), B(minor axis of the ellipse). Also plots the value of the mean.
        
    """


    def general_histograms(self, fcat):
        sns.set(style="white", palette="muted", color_codes=True)
        
        weights_mag_iso = np.ones_like(fcat['mag_iso'])/len(fcat['mag_iso'])
        weights_ellip = np.ones_like(fcat['ellipticity'])/len(fcat['ellipticity'])
        weights_a = np.ones_like(fcat['A'])/len(fcat['A'])
        weights_b = np.ones_like(fcat['B'])/len(fcat['B'])
        weights_flux_iso = np.ones_like(fcat['flux_iso'])/len(fcat['flux_iso'])
        
        sns.distplot(fcat['mag_iso'], axlabel = 'mag_iso', hist=True, hist_kws={'weights': weights_mag_iso})
        plt.axvline(fcat['mag_iso'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('$m(iso)$')
        plt.ylabel('$Frecuency$')
        plt.xlim(10, 30)
        plt.show()
        
        sns.distplot(fcat['ellipticity'], color="r", axlabel = 'ellipticity', hist=True, hist_kws={'weights': weights_ellip})
        plt.axvline(fcat['ellipticity'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('$ellipticity$')
        plt.ylabel('$Frecuency$')
        plt.show()
        
        sns.distplot(fcat['A'], color="y", bins = 20, axlabel = 'a', hist=True, hist_kws={'weights': weights_a})
        plt.axvline(fcat['A'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('$<A>$')
        plt.xlim(0, 50)
        plt.ylabel('$Frecuency$')
        plt.show()
        
        sns.distplot(fcat['B'], color="m", bins = 20,  axlabel = 'b', hist=True, hist_kws={'weights': weights_b})
        plt.axvline(fcat['B'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('$<B>$')
        plt.xlim(0, 10)
        plt.ylabel('$Frecuency$')
        plt.show()
    
    
        sns.distplot(fcat['flux_iso'], color="b", axlabel = 'flux_iso', hist=True, hist_kws={'weights': weights_flux_iso})
        plt.axvline(fcat['flux_iso'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('$F(iso)$')
        plt.ylabel('$Frecuency$')
        plt.show()


    """
        Define the elliptical-2DGaussian model to plot new celestial objects.
    
    """

    def gaussian2D(self, x, y, x_0, y_0, a, b, M):
        return M * math.exp(-0.5*((x-x_0)*(x-x_0)/(a*a) + (y-y_0)*(y-y_0)/(b*b)))
    
    
    """
        Define the relation between m(F) according to the astrophysical definition
        
    """
    
    def f_way_back(self,f,a=0, b=0):
        return -2.5*(math.log10(f/a) - b)
    
    """
        Transform the float-value of the x_position and y_position for celestial objects determined by Source Extractor into a integer value due to the fact that we requires to obtain a masking matrix which is discrete.
        
    """

    def transform_to_int(self, fcat):
        x_position = fcat['x']
        y_position = fcat['y']
        self.x_position_int= x_position.astype(int)   #transform double into int
        self.y_position_int= y_position.astype(int)
    
    """
        Method that opens and reads a .fits picture using Astropy library whose classes and methods are defined at http://astropy.readthedocs.org/en/latest/io/fits/index.html?highlight=fits#module-astropy.io.fits . Gets the value of the x_axis and y_axis as well as the intensity value at each pixel.
        
    """


    def open_read_picture(self, picture):
        hdulist_data_image=fits.open(picture, memmap=True)
        self.x_data_image=hdulist_data_image[0].header['NAXIS1']
        self.y_data_image=hdulist_data_image[0].header['NAXIS2']
        print 'The picture has a size of ({}x{})\n'.format(self.x_data_image, self.y_data_image)
        self.picture_data = hdulist_data_image[0].data
    
    """
        Method that creates a masked matrix only with 0 and 1. Attach the value 0 at those pixel where no celestial object was found by Source Extractor, and attach the value 1 where celestial objects where placed by Source Extractor.
        
    """

    def masking_matrix(self, picture):
        self.transform_to_int(self.fcat)
        self.open_read_picture(picture)
        self.matrix_data = np.zeros(shape=(self.y_data_image, self.x_data_image))
        for i in range (0, len(self.x_position_int)):
            self.matrix_data[self.y_position_int[i], self.x_position_int[i]]=1


    """
        Method that determined the percentage of 1 and 0 given at the masked matrix
    
    """


    def counts_percentages(self):
        self.total = self.x_data_image * self.y_data_image
        self.cont_1 = np.sum(self.matrix_data)
        self.cont_0 = self.total - self.cont_1

        self.percentage_1=(self.cont_1/self.total)*100.0
        self.percentage_0=(self.cont_0/self.total)*100.0
        print 'The total amount of pixels is :{}\n'.format(self.total)
        print 'The count of 1 is: {} \nThe count of 0 is: {} \nThe percentage of 1 is: {} \nThe percentage of 0 is: {}\n'.format(self.cont_1, self.cont_0, self.percentage_1, self.percentage_0)


    """
        Method that gets the number of new objects that we need to add to the picture if we want to get a percentage of 1 equals to a packing fraction called eta.
    
    """
    def packing_percentage(self, number_objects = 10):
        self.counts_percentages()
        #number_to_packing_double = eta * self.cont_1 /self.percentage_1
        #self.number_to_packing = int(number_to_packing_double)
        #self.packing = eta
        self.number_to_packing = number_objects
    
    
    """
        Method that writes and generates the new .fits picture that contains the new celestials objects using Astropy library at http://astropy.readthedocs.org/en/latest/io/fits/index.html?highlight=fits#module-astropy.io.fits
        
    """
    def get_simulation_picture(self, ID_number, PICTURE):
        self.simulation_picture = '{}_simulation_{}.fits'.format(PICTURE, ID_number)
        fits.writeto(self.simulation_picture, self.matrix_data)
    
    
    """
        Method that plots the elliptical celestial objects where no previous object was found (at those 0 matrix_data elements) randomly for a mag_iso value passed as attribute. Pick the value of A, B and ellipticity at the mean obtained from histograms. The relation between the intensity of the new picture pixels and the magnitude is made thanks to the exponential fitting.
        
    """
    
    def objectcreator_magnitude(self, mag_value = 0.0, n = 5, pic = ''):
        print 'Show the relation between the intensity and the manitude \n'
        
        #self.plot_fitting_exp(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        model, self.parameter_a, self.parameter_b = self.plot_fitting_exp_owndefined(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        
        intensity_value = model.f(m = mag_value, b = self.parameter_b, a = self.parameter_a)
        intensity_value_mag = model.f(m = self.mean_mag, b = self.parameter_b, a = self.parameter_a)
        
        print 'Value parameter_a is {}'.format(self.parameter_a)
        print 'Value parameter_b is {}'.format(self.parameter_b)
        print 'Value intensity is {}'.format(intensity_value)
        print 'mean_b is {}\nmean_a is {}'.format(self.mean_b, self.mean_a)
        print 'error_mean_b is {}\n error_mean_a is {}'.format(self.error_mean_b, self.error_mean_a)
        
        y_pixel = int(round(n*self.mean_b))
        x_pixel = int(round(n*self.mean_a))
        
        #print y_pixel
        #print x_pixel
        
        self.x_position_simulation = np.zeros(self.number_to_packing)
        self.y_position_simulation = np.zeros(self.number_to_packing)
        
        cont_percentage = 0
        
        while cont_percentage != self.number_to_packing:
            
            x = int(self.x_data_image * random.random())
            y = int(self.y_data_image * random.random())
            
            while x+x_pixel>self.x_data_image or x-x_pixel<0:
                x = int(self.x_data_image * random.random())
            
            while y+y_pixel>self.y_data_image or y-y_pixel<0:
                y = int(self.y_data_image * random.random())
            
            
            if self.matrix_data[y,x] == 0:
                
                #print 'value y {}'.format(y)
                #print 'value x {}'.format(x)
                
                self.x_position_simulation[cont_percentage] = x
                self.y_position_simulation[cont_percentage] = y
                
                self.matrix_data[y,x] = intensity_value
                
                for k in range (y-y_pixel, y+y_pixel):
                    for i in range (x-x_pixel, x+x_pixel):
                        if self.matrix_data[k,i]==0:
                            self.matrix_data[k,i]= self.picture_data[k,i]+ self.gaussian2D(i, k, x, y, self.mean_a, self.mean_b, intensity_value)
        
            cont_percentage = cont_percentage + 1

        #Attach intensity value to the rest of the picture pixels. In case you want to obtain only a pic with the new celestial objects please comment these lines
        for g in range (self.y_data_image):
            for f in range (self.x_data_image):
                if self.matrix_data[g,f]<=0.2:
                    self.matrix_data[g,f]=self.picture_data[g,f]
                elif self.matrix_data[g,f] == 1:
                    self.matrix_data[g,f]=self.picture_data[g,f]

        #write in a new picture the matrix_data
        self.get_simulation_picture(mag_value, pic)



    """
        Method that looks for the celestial objects already created by the previous method in the new catalog obtained after running Source Extractor.
        
        PROBLEM: super inefficient due to the fact of using Python for-loops
    
    """

    def searcher(self, fcat, fcat_simulation):
        
        print '\nTry searcher\n'
        
        x_position_fcat_simulation = fcat_simulation['x'].astype(int)
        y_position_fcat_simulation = fcat_simulation['y'].astype(int)
        
        self.out_mag = np.zeros(self.number_to_packing) #define out_mag
        
        #print len(self.x_position_simulation)
        #print len(x_position_fcat_simulation)
        
        length_fcat_simulation = len(x_position_fcat_simulation)
        length_x_position_simulation = len(self.x_position_simulation)
        
        for i in range (0, length_fcat_simulation):
            for k in range (0, length_x_position_simulation):
                if (x_position_fcat_simulation[i] == self.x_position_simulation[k] and y_position_fcat_simulation[i] == self.y_position_simulation[k]):
                    self.out_mag[k] = fcat_simulation['mag_iso'][i]
                    print self.out_mag[k]
                    break
        print self.out.mag




    """
        Method that looks for the celestial objects already created by the previous method in the new catalog obtained after running Source Extractor.
        
        TRY: diccionaries with a double key-tag
        
        PROBLEM: quite fast, only int-keys
        
        OPTIONS in FUTURE: may it be implementated in C++ using mapping
    
    """
    
    
    
    def searcher_dic(self, fcat, fcat_simulation):
        
        print '\nTry dic\n'
        
        #Define super dic empty
        d={}
        #Get values of x and y from new catag
        x_position_fcat_simulation = fcat_simulation['x'].astype(int)
        y_position_fcat_simulation = fcat_simulation['y'].astype(int)
        
        #self.out_mag_all = np.zeros(self.number_to_packing) #define out_mag
        
        #print len(self.x_position_simulation)
        #print len(x_position_fcat_simulation)
        
        length_fcat_simulation = len(x_position_fcat_simulation)
        length_x_position_simulation = len(self.x_position_simulation)
        
        for i in range (0, length_fcat_simulation):
            #print type(x_position_fcat_simulation[i])
            d[(int(x_position_fcat_simulation[i]), int(y_position_fcat_simulation[i]))] = fcat_simulation['flux_iso'][i]
        
        for k in range (0, length_x_position_simulation):
            flux = d.get((int(self.x_position_simulation[k]), int(self.y_position_simulation[k])), -1)
            if flux is not -1:
                self.out_mag.append(self.f_way_back(flux,a=self.parameter_a,b=self.parameter_b))


    """
        Method that looks for the celestial objects already created by the previous method in the new catalog obtained after running Source Extractor.
        
        TRY: KDTree altorigthm
        
        ADVANTAGES: really fast
    
    """
        
    def searcher_kdtree(self, fcat, fcat_simulation, FILE_NAME):
        
        print '\nTry KDTree'
        x_position_fcat_simulation = fcat_simulation['x']
        y_position_fcat_simulation = fcat_simulation['y']
        position_all = zip(x_position_fcat_simulation.ravel(), y_position_fcat_simulation.ravel())
        tree = spatial.KDTree(position_all)
        position_obj_created = zip(self.x_position_simulation, self.y_position_simulation)
                        
        self.posible_obj_distances, self.posible_obj_index = tree.query(position_obj_created, distance_upper_bound = 3*self.mean_a)

        cont_lost_obj = 0

        for index in self.posible_obj_index:
            if index < len(position_all):
                self.out_mag.append(fcat_simulation['mag_iso'][index])
                self.out_flux.append(fcat_simulation['flux_iso'][index])
                self.out_flux_max.append(fcat_simulation['flux_max'][index])
            #self.out_mag_after_transf.append(self.f_way_back(fcat_simulation['flux_iso'][index],a=112.188243402,b=9.95005045126))
            if index == len(position_all):
                cont_lost_obj = cont_lost_obj + 1
    
        self.lost_objects.append(cont_lost_obj)

        print 'The number of lost galaxies is {}\n'.format(cont_lost_obj)

