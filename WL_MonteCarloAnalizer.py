# Name: WL_MonteCarlo.py
#
# TFG Program VII
#
# Description: Prototype of the script requires to deal with the noise
#
# Returns: fits image corresponding to the final mass-density map

__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "2.0.0"
__maintainer__ = "Guadalupe Canas"
__email__ = "gch24@alumnos.unican.es"


#
# Improvements: transformed to be callable function from WL_script.py or CatalogPlotter3.py. Also included how to make all cuts at the same time
#

import matplotlib.pyplot as plt #Works for making python behave like matlab
import numpy as np #Maths arrays and more
from scipy.special import erfcinv #Inverse Complementary Error Function
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
from astropy.io import fits #Open and Reading Files
import matplotlib.pylab as P #For histograms
import matplotlib
import seaborn as sns
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

#
# Brief description of the program: we need to read 1000 randomized images x number of cuts. Example: 1000 randomized images of a cut from 12 until 18 mag. If we have 4 cuts, then we have 4000 randomized images to read.
#


def sigma_maker(data_image, cycles, local):
    matplotlib.use('GTK')
    #data_image is the output image that stay always the same from
    hdulist_data_image=fits.open(data_image, memmap=False)
    #(1): Required to know how many pixels both images (are the same values so we can just work with one). These information is saved in NAXIS1 and NAXIS2
    x_data_image=hdulist_data_image[0].header['NAXIS1']
    y_data_image=hdulist_data_image[0].header['NAXIS2']
    
    #(2): Obtaining data from FITS images
    sci_data_image = hdulist_data_image[0].data
    
    #Need to filter that array to obtain only one layer
    sci_data_image_st=sci_data_image[0,:,:]
    
    #(3): Prepare for -for loop- setting matrix into zero
    matrix=np.zeros(shape=(y_data_image, x_data_image))
    S_N=np.zeros(shape=(y_data_image, x_data_image))
    #print matrix.shape
    matrix_sum=np.zeros(shape=(y_data_image, x_data_image))
    hdulist_data_image.close(closed=True)
    list_hdu=[''] * cycles #list for saving data of randomized images
    inv_b=1.0/cycles
    #print inv_b
    histogram_array_0=np.zeros(cycles)
    for k in range (0, 81):
        for j in range(0, 100):
            for i in range (0, cycles):
                print i
                #std_{}_{} is the output randomized image that changes
                hdulist_random_image=fits.open('std_{}_{}.fits'.format(local, i), memmap=False)
                list_hdu[i]=hdulist_random_image[0].data
            
                histogram_array_0[i]=list_hdu[i].item((k,j))
                hdulist_random_image.close(closed=True)
    
                mu, std = norm.fit(histogram_array_0)
                matrix[k,j]=std
    
                print std
    
    S_N = np.divide(sci_data_image_st, matrix)
  
    S_N[S_N > 4.0] = S_N[S_N > 4.0]
    
    #(9): We plot this information in the new file
    fits.writeto('WL_map4sigma_{}.fits'.format(local), S_N)

