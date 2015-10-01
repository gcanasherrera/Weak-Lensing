# Name: Class_ObjectCreator.py
#
# Bachelor Disertation Program VI
#
# Type: python class
#
# Content: 1 Class, 1 constructor,
#
# Description: General Class destinated to produce a simulation over a .fits picture in order to determinate whether Source Extractor
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "1.0.0"
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
import random #pseudo-random generator
import seaborn as sns #Improvements for statistical-plots
from pymodelfit import FunctionModel1DAuto #Create own model of fitting


class magnitude_exponential(FunctionModel1DAuto):
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
        self.mean_intensity = 0.0
        self.mean_ellip = np.mean(fcat['ellipticity'])
        self.mean_b = np.mean(fcat['B'])
        self.mean_ec = np.sqrt(1-(1-self.mean_ellip)*(1-self.mean_ellip))
        self.mean_a = self.mean_b/math.sqrt(1-self.mean_ec*self.mean_ec)


        self.out_mag_all = []
        self.out_mag = []


        self.random_mag = []
        self.random_intensity = []
        self.random_ellip = []
        self.random_b = []
        self.random_a = []
    
    
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
        
        #plt.errorbar(x, y, self.fcat['fluxerr_iso'], self.fcat['magger_iso'], fmt='ko', label='Source Extractor')
        #plt.plot(x, y, 'ko', label='Source Extractor')
        xp= np.linspace(8.2, 25, len(x))
        plt.loglog(xp, t(xp), 'r--', label='Exponential Fitting')
        plt.xlabel('Mag_iso')
        plt.ylabel('Flux_iso')
        plt.legend()
        plt.title('Mag_iso Vs. Flux_iso')
        plt.show()
        print 'The fit of Flux_iso Vs. mag_iso is: {}\n'.format(t)
        #print t.alpha
        self.mean_intensity = t(self.mean_mag)
        print 'The mean intensity is : {}\nThe mean isophotal magnitude is : {} '.format(self.mean_intensity, self.mean_mag)
    
    
    
    """
        Method that fits the magnitude as a function of the flux in normal scale through an exponential model defined myself using
        
        """
    
    def plot_fitting_exp_owndefined(self, x, y):
        sns.set(style="white", palette="muted", color_codes=True)
        model = magnitude_exponential()
        af,bf = model.fitData(x,y)
        model.plot( logplot='xy')
        plt.text(0,11,'My Model:',fontsize=16)
        for i,k in enumerate(model.pardict):
            plt.text(0,10-i,'%s = %f'%(k,model.pardict[k]),fontsize=16)
        plt.xlabel('Mag_iso')
        plt.ylabel('Flux_iso')
        
        plt.show()
        self.mean_intensity = model.f(m = self.mean_mag, a = af, b = bf)
        print 'The mean intensity is : {}\nThe mean isophotal magnitude is : {} '.format(self.mean_intensity, self.mean_mag)
        
        
        return model, af, bf
    
    
    
    
    
    """
        Method that plots the histograms corresponding to the magnitudes: flux_iso, mag_iso, ellipticity, A (major axis of the ellipse), B(minor axis of the ellipse). Also plots the value of the mean.
        
    """


    def general_histograms(self, fcat):
        sns.set(style="white", palette="muted", color_codes=True)
        rs = np.random.RandomState(10)
        
        sns.distplot(fcat['flux_iso'], color="b", axlabel = 'flux_iso', hist=True)
        plt.axvline(fcat['flux_iso'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.show()
        
        sns.distplot(fcat['mag_iso'], axlabel = 'mag_iso', hist=True)
        plt.axvline(fcat['mag_iso'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.show()
        
        sns.distplot(fcat['ellipticity'], color="r", axlabel = 'ellipticity', hist=True)
        plt.axvline(fcat['ellipticity'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.show()
        
        sns.distplot(fcat['A'], color="y", axlabel = 'a', hist=True)
        plt.axvline(fcat['A'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.show()
        
        sns.distplot(fcat['B'], color="m", axlabel = 'b', hist=True)
        plt.axvline(fcat['B'].mean(), color='k', linestyle='dashed', linewidth=2)
        plt.show()


    """
        Define the elliptical-2DGaussian model to plot new celestial objects.
    
    """

    def gaussian2D(self, x, y, x_0, y_0, a, b, M):
        return M * math.exp(-0.5*((x-x_0)*(x-x_0)/(a*a) + (y-y_0)*(y-y_0)/(b*b)))
    
    
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
        hdulist_data_image=fits.open(picture, memmap=False)
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
    def packing_percentage(self, eta = 0.08):
        self.counts_percentages()
        number_to_packing_double = eta * self.cont_1 /self.percentage_1
        self.number_to_packing = int(number_to_packing_double)
        self.packing = eta
    
    
    """
        Method that writes and generates the new .fits picture that contains the new celestials objects using Astropy library at http://astropy.readthedocs.org/en/latest/io/fits/index.html?highlight=fits#module-astropy.io.fits
        
    """
    def get_simulation_picture(self, ID_number):
        self.simulation_picture = '{}_Simulation_{}.fits'.format('w2_53_stack', ID_number)
        fits.writeto(self.simulation_picture, self.matrix_data)
    
    
    """
        Method that plots the elliptical celestial objects where no previous object was found (at those 0 matrix_data elements) randomly for a mag_iso value passed as attribute. Pick the value of A, B and ellipticity at the mean obtained from histograms. The relation between the intensity of the new picture pixels and the magnitude is made thanks to the exponential fitting.
        
    """
    
    def objectcreator_magnitude(self, mag_value = 0.0, n = 5):
        print 'Show the relation between the intensity and the manitude \n'
        
        self.plot_fitting_exp(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        model, af, bf = self.plot_fitting_exp_owndefined(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        #intensity_value = fitting_intensity_mag(mag_value)
        
        intensity_value = model.f(m = mag_value, b = bf, a = af)
        intensity_value_mag = model.f(m = self.mean_mag, b = bf, a = af)
        
        print af
        print bf
        
        print 'Value intensity is {}'.format(intensity_value)
        print 'mean value intensity is {}'.format(intensity_value_mag)
        
        print 'mean_b is {}\nmean_a is {}'.format(self.mean_b, self.mean_a)
        
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
                
                for k in range (y-y_pixel, y+y_pixel):
                    for i in range (x-x_pixel, x+x_pixel):
                        if self.matrix_data[k,i]==0:
                            self.matrix_data[k,i]= self.gaussian2D(i, k, x, y, self.mean_a, self.mean_b, intensity_value)
        
            cont_percentage = cont_percentage + 1

        #Attach intensity value to the rest of the picture pixels.
        for g in range (self.y_data_image):
            for f in range (self.x_data_image):
                if self.matrix_data[g,f]<=0.1:
                    self.matrix_data[g,f]=self.picture_data[g,f]

        #write in a new picture the matrix_data
        self.get_simulation_picture(mag_value)








    """
        Method that looks for the celestial objects already created by the previous method in the new catalog obtained after running Source Extractor.
        
        PROBLEM: super inefficient due to the fact of using Python for-loops
    
    """

    def searcher(self, fcat, fcat_simulation):
        
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
            d[(int(x_position_fcat_simulation[i]), int(y_position_fcat_simulation[i]))] = fcat_simulation['mag_iso'][i]
        
        #print d
        
        for k in range (0, length_x_position_simulation):
            mag = d.get((int(self.x_position_simulation[k]), int(self.y_position_simulation[k])), -1)
            if mag is not -1:
                self.out_mag.append(mag)
            
#print self.out_mag



    """
    Method that plots the elliptical celestial objects where no previous object was found (at those 0 matrix_data elements) randomly for a ellipticity value passed as attribute. Pick the value of A, B and ellipticity at the mean obtained from histograms. The relation between the intensity of the new picture pixels and the magnitude is made thanks to the exponential fitting.
    
    """
        
    def objectcreator_ellipticity(self, ellip_value = 0.0, n = 5):
        print 'Show the relation between the intensity and the manitude \n'
        
        fitting_intensity_mag = self.plot_fitting_exp(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        intensity_value = fitting_intensity_mag(self.mean_mag)
        
        print 'mean_b is {}\nmean_a is {}'.format(self.mean_b, self.mean_a)
        
        y_pixel = int(round(n*self.mean_b))
        x_pixel = int(round(n*self.mean_a))
        
        #print y_pixel
        #print x_pixel
        
        self.x_position_simulation = np.zeros(self.number_to_packing)
        self.y_position_simulation = np.zeros(self.number_to_packing)
        
        cont_percentage = 0
        
        
        self.mean_ec = np.sqrt(1-(1-self.ellip_value)*(1-self.ellip_value))
        self.mean_a = self.mean_b/math.sqrt(1-self.mean_ec*self.mean_ec)
        
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
                
                for k in range (y-y_pixel, y+y_pixel):
                    for i in range (x-x_pixel, x+x_pixel):
                        if self.matrix_data[k,i]==0:
                            self.matrix_data[k,i]= self.gaussian2D(i, k, x, y, self.mean_a, self.mean_b, intensity_value)

        cont_percentage = cont_percentage + 1

        #Attach intensity value to the rest of the picture pixels.
        for g in range (self.y_data_image):
            for f in range (self.x_data_image):
                if self.matrix_data[g,f]<=0.1:
                    self.matrix_data[g,f]=self.picture_data[g,f]
                
        #write in a new picture the matrix_data
        self.get_simulation_picture(mag_value)


