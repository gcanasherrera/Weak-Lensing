# Name: filter_mag_gal.py
#
# Bachelor Disertation Program IV
#
# Description: Python File containing a function that splits the catalog of galaxies (already corrected) in different catalog cuts, creating a mass-density map FITS file (data_image and random_image) given a number of files.
#
# For further information about the functions here defined the users are kindly referred to the DOCUMENTATION file attached in this package.
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "1.1.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"


# Improvements: Part of CatalogPlotter.3 ---> include

import matplotlib.pyplot as plt #Works for making python behave like matlab
#import sextutils as sex #Program used to read the original catalog
import numpy as np #Maths arrays and more
import numpy.ma as ma #Masking arrays
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
#from astropy.modeling import models, fitting #Package for fitting Legendre Polynomials
import warnings #Advices
from mpl_toolkits.mplot3d import Axes3D #Plotting in 3D
#import ellip_fitting3 as ellip_fit #Ellipticity fitting
#from CatalogPlotter3 import ellipticity
#from CatalogPlotter3 import plotter


def filter_mag(catalog, init, fin):
    BEFORE_NAME = catalog.find('.')
    FILE_NAME = catalog[:BEFORE_NAME]

    names_ellipto = ["x", "y", "mag_iso", "median", "ixx", "iyy", "ixy", "a_input", "b_input", "theta", "ellipticity", "errcode", "sigsky", "size", "flux", "mean_rho_4th", "sigma_e", "wander"]

    print("")
    print("Now we make cuts of the number of galaxies as a function of the magnitudes")
    print("")
    print("Let's obtain cuts of the galaxies catalog as a function of mag. We need to bound")
    print("")
    fcat_shapes_galaxies_corrected= np.genfromtxt(catalog, names=names_ellipto)
    #Look for maximum and minimum of mag_iso
    mag_iso_max_div=np.amax(fcat_shapes_galaxies_corrected['mag_iso'])
    mag_iso_min_div=np.amin(fcat_shapes_galaxies_corrected['mag_iso'])

    #plotter(fcat_shapes_galaxies_corrected, "mag_iso", "size", 1)
    plt.show(block=False)
    print('The minimum of mag_iso is {} and the maximum of mag_iso is {}'.format(mag_iso_min_div, mag_iso_max_div))
    #groups=input("How many groups you want to have?: ", 3)
    groups=3
    #cut_values=np.zeros(groups-1)
    #for k in range (0, groups-1):
     #  print('Enter the value of the cut {}: '.format(k))
      # cut=input()
       #cut_values[k]=cut
    
    cut_values=np.array([20, 22])
    print 'The introduced values for the cuts are: {} and the length {}'.format(cut_values, len(cut_values))
    # First cut taking into account the minimum

    cont=0
    k=0
    
    for i in range (init, fin):
        
        print('\n Round : {} \n'.format(i))
        
        catalog_name_mask_min= 'g_min_{}.fcat'.format(FILE_NAME)
        terminal_galaxies_mask_min= 'perl fiatfilter.pl -v "mag>{} && mag<{}" {}>{}'.format(mag_iso_min_div, cut_values[0], catalog, catalog_name_mask_min)
        subprocess.call(terminal_galaxies_mask_min, shell=True)
        
        #print catalog_name_mask_min
        terminal_fiatmap_min='./fiatmap -b 25 -w /gpfs/csic_users/canasg/area  -r std_min_{:n}.fits {} 300 3000 g_min_{}_{:n}.fits'.format(i, catalog_name_mask_min, FILE_NAME, i)
        subprocess.call(terminal_fiatmap_min, shell=True)
        # Last cut taking into account the maximum
        catalog_name_mask_max= 'g_max_{}.fcat'.format(FILE_NAME)
        terminal_galaxies_mask_max= 'perl fiatfilter.pl -v "mag>{} && mag<{}" {}>{}'.format(np.amax(cut_values), mag_iso_max_div, catalog, catalog_name_mask_max)
        
        #Read values of the new mask catalogs
        #fcat_mask_min = np.genfromtxt(catalog_name_mask_min, names=names_ellipto)
        #Call the function ellipticity
        #ellipticity(fcat_mask_min, 7)
        #plt.show()
        
        #Now we fix in case that the length of cut_values is different than one
        if len(cut_values)!=1:
            for a in range (0, len(cut_values)):
                if a+1<len(cut_values):
                    catalog_name_mask='g_{}_{}.fcat'.format(a, FILE_NAME)
                    terminal_galaxies_mask= 'perl fiatfilter.pl -v "mag>{} && mag<{}" {}>{}'.format(cut_values[a], cut_values[a+1], catalog, catalog_name_mask)
                    subprocess.call(terminal_galaxies_mask, shell=True)
                    terminal_fiatmap='./fiatmap -b 25  -w  /gpfs/csic_users/canasg/area -r std_{:n}_{:n}.fits {} 300 3000 g_{:n}_{}_{:n}.fits'.format(a, i, catalog_name_mask, a, FILE_NAME, i)
                    subprocess.call(terminal_fiatmap, shell=True)
                    #Read values of the new mask catalogs
                    fcat_mask = np.genfromtxt(catalog_name_mask, names=names_ellipto)
                    #plt.figure()
                    #ellipticity(fcat_mask, 8+a)
                    plt.show(block=False)

        subprocess.call(terminal_galaxies_mask_max, shell=True)
        fcat_mask_max = np.genfromtxt(catalog_name_mask_max, names=names_ellipto)
        terminal_fiatmap_max='./fiatmap -b 25 -w /gpfs/csic_users/canasg/area -r std_max_{:n}.fits {} 300 3000 g_max_{}_{:n}.fits'.format(i, catalog_name_mask_max, FILE_NAME, i)
        subprocess.call(terminal_fiatmap_max, shell=True)
        #ellipticity(fcat_mask_max,14)
        plt.show()
        cont=cont + 1

    print("END !!!")
    print("")
    return groups
