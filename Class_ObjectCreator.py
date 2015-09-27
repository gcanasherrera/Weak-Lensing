# Name: ObjectCreator.py
#
# Bachelor Disertation Program VI
#
# Type: Class
#
# Description:
#
# Returns: FITS image including the objects we created
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
from astropy.io import fits #Open and Reading Files
import random #pseudo-random generator
import seaborn as sns #Improvements for statistical-plots

class ObjectCreator(object):

    """
    Creates c
    
    """
    #define attributes of future class-objects
    
    def __init__(self, fcat):
        
        self.fcat = fcat
        self.simulation_picture = ''
    
        self.mean_mag = np.mean(fcat['mag_iso'])
        self.mean_intensity = 0.0
        self.mean_ellip = np.mean(fcat['ellipticity'])
        self.mean_b = np.mean(fcat['B'])
        self.mean_ec = np.sqrt(1-(1-self.mean_ellip)*(1-self.mean_ellip))
        self.mean_a = self.mean_b/math.sqrt(1-self.mean_ec*self.mean_ec)


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

    def plot_fitting_exp(self, x, y):
        sns.set(style="white", palette="muted", color_codes=True)
        t_init = models.ExponentialCutoffPowerLaw1D(amplitude=10., x_0=10., alpha=100., x_cutoff=10)
        fit_t = fitting.LevMarLSQFitter()
        t = fit_t(t_init, x, y)
        plt.figure()
        plt.plot(x, y, 'ko', label='Source Extractor')
        xp= np.linspace(8.2, 25, len(x))
        plt.plot(xp, t(xp), 'r--', label='Exponential Fitting')
        plt.xlabel('Mag_iso')
        plt.ylabel('Flux_iso')
        plt.legend()
        plt.title('Mag_iso Vs. Flux_iso')
        plt.show()
        print 'The fit of Flux_iso Vs. mag_iso is: {}\n'.format(t)
        print t.alpha
        self.mean_intensity = t(self.mean_mag)
        print 'The mean intensity is : {}\nThe mean isophotal magnitude is : {} '.format(self.mean_intensity, self.mean_mag)
        
        return t


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


    def gaussian2D(self, x, y, x_0, y_0, a, b, M):
        return M * math.exp(-0.5*((x-x_0)*(x-x_0)/(a*a) + (y-y_0)*(y-y_0)/(b*b)))
    

    def transform_to_int(self, fcat):
        x_position = fcat['x']
        y_position = fcat['y']
        self.x_position_int= x_position.astype(int)   #transform double into int
        self.y_position_int= y_position.astype(int)
    

    def open_read_picture(self, picture):
        hdulist_data_image=fits.open(picture, memmap=False)
        self.x_data_image=hdulist_data_image[0].header['NAXIS1']
        self.y_data_image=hdulist_data_image[0].header['NAXIS2']
        print 'The picture has a size of ({}x{})\n'.format(self.x_data_image, self.y_data_image)
        self.picture_data = hdulist_data_image[0].data

    def masking_matrix(self, picture):
        self.transform_to_int(self.fcat)
        self.open_read_picture(picture)
        self.matrix_data = np.zeros(shape=(self.y_data_image, self.x_data_image))
        for i in range (0, len(self.x_position_int)):
            self.matrix_data[self.y_position_int[i], self.x_position_int[i]]=1

    def counts_percentages(self):
        self.total = self.x_data_image * self.y_data_image
        self.cont_1 = np.sum(self.matrix_data)
        self.cont_0 = self.total - self.cont_1

        self.percentage_1=(self.cont_1/self.total)*100.0
        self.percentage_0=(self.cont_0/self.total)*100.0
        print 'The total amount of pixels is :{}\n'.format(self.total)
        print 'The count of 1 is: {} \nThe count of 0 is: {} \nThe percentage of 1 is: {} \nThe percentage of 0 is: {}\n'.format(self.cont_1, self.cont_0, self.percentage_1, self.percentage_0)

    def packing_percentage(self, eta = 0.08):
        self.counts_percentages()
        number_to_packing_double = eta * self.cont_1 /self.percentage_1
        self.number_to_packing = int(number_to_packing_double)
        self.packing = eta
    
    
    def get_simulation_picture(self, ID_number):
        self.simulation_picture = '{}_Simulation_{}.fits'.format('w2_53_stack', ID_number)
        fits.writeto(self.simulation_picture, self.matrix_data)
    
    
    def objectcreator_mean(self, mag_value = 0.0, n = 5):
        print 'Show the relation between the intensity and the manitude \n'
        
        fitting_intensity_mag = self.plot_fitting_exp(self.fcat['mag_iso'], self.fcat['flux_iso'])
        
        intensity_value = fitting_intensity_mag(mag_value)
        
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







