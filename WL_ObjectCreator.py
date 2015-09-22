# Name: ObjectCreator.py
#
# Bachelor Disertation Program VI
#
# Description:
#
# Returns: FITS image including the objects we created
#

__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "1.0.0"
__maintainer__ = "Guadalupe Canas"
__email__ = "gch24@alumnos.unican.es"

import matplotlib.pyplot as plt #Works for making python behave like matlab
#import sextutils as sex #Program used to read the original catalog
import numpy as np #Maths arrays and more
import numpy.ma as ma #Masking arrays
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
from astropy.modeling import models, fitting #Package for fitting Legendre Polynomials
import warnings #Advices
from astropy.io import fits #Open and Reading Files
import random
import matplotlib.pylab as P #For histograms
import matplotlib
import seaborn as sns
from scipy.optimize import curve_fit

#Define elliptical-gaussial2D

def gaussian2D(x, y, x_0, y_0, a, b, M):
    return M * math.exp(-0.5*((x-x_0)*(x-x_0)/(a*a) + (y-y_0)*(y-y_0)/(b*b)))

def plot_fitting_pol(mag_1, mag_2, val_1, val_2, lin):
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
    plt.plot(np.log(mag_1), np.log(mag_2), 'ro', label='Source Extractor')
    #plt.errorbar(np.log(fcat['flux_iso']), np.log(fcat['mag_iso']), xerr = xerr_norm, yerr = yerr_norm, fmt='r.')
    plt.plot(xp, p3(xp), 'k-', label='fitting')
    plt.legend()
    plt.show(block=False)
    return p3

def plot_fitting_exp(x, y):
    t_init = models.ExponentialCutoffPowerLaw1D(amplitude=10., x_0=10., alpha=100., x_cutoff=10)
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, x, y)
    plt.figure()
    plt.plot(x, y, 'ko', label='Source Extractor')
    xp= np.linspace(10, 25, len(x))
    plt.plot(x, t(x), 'r.', label='fitting')
    plt.xlabel('Mag_iso')
    plt.ylabel('Flux_iso')
    plt.legend()
    plt.title('Mag_iso Vs. Flux_iso')
    plt.show(block=False)
    print 'The fit of Flux_iso Vs. mag_iso is: {}\n'.format(t)
    return t


def main():
    #We need to create several objects with random position, magnitude, ellipticity and Orientation

    print("ObjectCreator introduces celestial objects in the fits image")
    print("")
    picture= 'w2_53_stack.fits'
    #picture = raw_input("Please, introduce the name of the fits image you want to read: ")
    #catalog = raw_input("Please, introduce the name of catalog obtained from Source Extractor: ")
    catalog='w2_53_stack.fcat'
    catalog_fiat = 'fiat_{}'.format(catalog)
    transform_into_fiat='perl sex2fiat.pl {}>{}'.format(catalog, catalog_fiat)
    subprocess.call(transform_into_fiat, shell=True)
    names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
    fcat = np.genfromtxt(catalog_fiat, names=names)
        
    # (0.1): Fits the mag_iso againts the flux_iso
    p3_mag_vs_flux = plot_fitting_pol(fcat['flux_iso'], fcat['mag_iso'], 'flux_iso', 'mag_iso', 22)
    # (0.2); Fits flux_iso againts the  mag_iso
    fit_exp_flux_vs_mag = plot_fitting_exp(fcat['mag_iso'], fcat['flux_iso'])
    
    
    
    
    plt.figure()
    sns.set_style("white")
    sns.set_style("ticks")
    plt.title('Histogram a')
    plt.xlabel('a')
    plt.ylabel('Frequency')
    P.hist(fcat['A'], 50, normed=1, histtype='stepfilled')
    P.show(block=False)
    
    plt.figure()
    plt.title('Histogram Flux_Max')
    plt.xlabel('Flux_max')
    plt.ylabel('Frequency')
    P.hist(fcat['flux_max'], 50, normed=1, histtype='stepfilled')
    P.show(block=False)
    
    # (0.2): Plots histograms
    sns.set(style="white", palette="muted", color_codes=True)
    rs = np.random.RandomState(10)

    # Set up the matplotlib figure
    f, axes = plt.subplots(2, 2, figsize=(7, 7), sharex=False)
    #sns.despine(left=True)

    # Plot a simple histogram with binsize determined automatically
    sns.distplot(fcat['mag_iso'], kde=False, color="b", ax=axes[0, 0], axlabel = 'mag_iso')

    # Plot a kernel density estimate and rug plot
    sns.distplot(fcat['ellipticity'], kde=False, color="r", ax=axes[0, 1], axlabel = 'ellipticity')

    # Plot a filled kernel density estimate
    sns.distplot(fcat['A'], kde=False, color="g", ax=axes[1, 0], axlabel = 'a')

    # Plot a historgram and kernel density estimate
    sns.distplot(fcat['B'], kde=False, color="m", ax=axes[1, 1], axlabel = 'b')

    plt.setp(axes, yticks=[])
    plt.tight_layout()
    
    # (1): Determine position of the celestial object according to the catalog
    x_position = fcat['x']
    y_position = fcat['y']
    x_position_int= x_position.astype(int)   #transform double into int
    y_position_int= y_position.astype(int)
    plt.show()
    
    
    # (2): Create a mask of 1 (where there are objects) and 0 (where there are no objects). We need to pass a catalog obtained from Source Extractor whose threshold is quite low. --> La posicion de los objectos celestes tienen decimales y la matrix es discreta....

    hdulist_data_image=fits.open(picture, memmap=False)
    x_data_image=hdulist_data_image[0].header['NAXIS1']
    y_data_image=hdulist_data_image[0].header['NAXIS2']
    print 'The picture has a size of ({}x{})'.format(x_data_image, y_data_image)

    cont_1=0.0
    cont_0=0.0
    
    matrix_data=np.zeros(shape=(x_data_image, y_data_image))
    
    for i in range (0, len(x_position_int)):
        matrix_data[x_position_int[i], y_position_int[i]]=1

    # (3): Count the number of 1 and 0 and make percentages thanks to cont_1 and cont_2
    for j in range (0, x_data_image):
        for k in range (0, y_data_image):
            if matrix_data[j,k]==1:
                cont_1=cont_1+1.0
            else:
                cont_0=cont_0+1.0

    np.savetxt('test.txt', matrix_data, delimiter=',')

    total= x_data_image * y_data_image
    print total
    percentage_1=(cont_1/total)*100.0
    percentage_0=(cont_0/total)*100.0
    print 'The count of 1 is: {} \n The count of 0 is: {} \nThe percentage of 1 is: {} \nThe percentage of 0 is: {}'.format(cont_1, cont_0, percentage_1, percentage_0)
    
    # (4): Obtain the value of object we need to add in order to obtain a porcentage of 1 around 50%. We need to create an array with len = number_to_fifty
    
    number_to_fifty = 50 * cont_1 /percentage_1
    
    # (5): Make relationship between the value of intensity in the pixel and the value of the magnitude obtained from source extractor --> PROBLEMA: estableces esto en funcion de lo que te ha dado source extractor anteriormente... no tiene poco sentido?
    
    sci_data_image = hdulist_data_image[0].data #intensity data of the image
    intensity_relation = sci_data_image[x,y]
    mag_relation = 30.0
    mag_1_relation = intensity_relation/mag_relation
    
    # (5): Create arrays with random numbers for magnitude, ellipticity.
    random_mag = 30.0 * np.random.rand(number_to_fifty)
    random_mag_intensity = random_mag * mag_1_relation
    random_ellip = np.random.rand(number_to_fifty)
    
    # (6): Go from matrix_data to matrix_pixel
    matrix_pixel = np.zeros(shape=(x_data_image, y_data_image))
    picture_data = hdulist_data_image[0].data
    matrix_pixel[matrix_data == 1] = picture_data
    
    # (7): Loop for adding objects
    cont_fifty = 0 #contador to know how many object we have added
    a = b/math.sqrt(1-e*e)
    while cont_fifty != number_to_fifty:
        x = x_data_image * np.random.rand(0)
        y = y_data_image * np.random.rand(0)
        if matrix_pixel[x,y] == 0:
            matrix_pixel[x,y] = random_mag_intensity[cont_fifty]
            cont_fifty = cont_fifty + 1
            
            for i in range (x-5*a_mean, x+5*a_mean):
                for k in range (y-5*b_mean, x+5*b_mean):
                    matrix_pixel[i,k]=gaussian2D(i, k, x, y, a_mean, b_mean, M)
        
        else:
            continue

# COMO HACES QUE EL PIXEL NO SEA CUADRADO? PROPOSICIoN












if __name__ == "__main__":
    main()