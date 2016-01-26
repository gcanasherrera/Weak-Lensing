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
    #print matrix.shape
    matrix_sum=np.zeros(shape=(y_data_image, x_data_image))
    hdulist_data_image.close(closed=True)
    list_hdu=[''] * cycles #list for saving data of randomized images
    inv_b=1.0/cycles
    #print inv_b
    histogram_array_0=np.zeros(cycles)
    for i in range (0, cycles):
        print i
	#std_{}_{} is the output randomized image that changes
        hdulist_random_image=fits.open('std_{}_{}.fits'.format(local, i), memmap=False)
        list_hdu[i]=hdulist_random_image[0].data
        histogram_array_0[i]=list_hdu[i].item(5)
    	hdulist_random_image.close(closed=True)
	
        #(4): Ready for -masking- to go through the images at the same time. Set matrix_0=0
        matrix[sci_data_image_st < list_hdu[i]] = 1.0
        matrix_sum=matrix_sum + matrix
    
    P.figure()
    P.hist(histogram_array_0, 50, normed=1, histtype='stepfilled')
    P.show()
    
    print matrix_sum
    #(4): Need to create a new matrix_P whose elements are just P=elements(matrix_sum)/b
    matrix_P=np.multiply(matrix_sum, inv_b)

    #(5): Now we go from P into sigma using sigma=erfc^(-1)(2P)
    sigma=erfcinv(2*matrix_P)
    print("")
    
    print sigma
    
    #(6): Write back the data into a new .FITS image
    fits.writeto('WL_map_{}.fits'.format(local), sigma)
    
    #(7): we need to force our sigma matrix to be the size (306, 348) (this values are just set by hand)
    #rescale_sigma=np.zeros(shape=(306, 348))
    #rescale_sigma=sigma[9:88, 11:106]
    #print rescale_sigma.shape
 
    #(8): We need to analyze the matrix SIGMA in order to find the pixel whose values are > 4
    four_sigma=np.zeros(shape=(y_data_image, x_data_image))
    #print four_sigma.shape
    
    #four_sigma[rescale_sigma > 4.0] = rescale_sigma[rescale_sigma > 4.0]
    four_sigma[sigma > 4.0] = sigma[sigma > 4.0]
    
    #(9): We plot this information in the new file
    fits.writeto('WL_map4sigma_{}.fits'.format(local), four_sigma)

    return

