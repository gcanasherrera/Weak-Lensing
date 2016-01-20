# Name: WL_FitsSubstractor.py
#
# TFG Program V
#
# Description: reads two fits files and substracts intensity between them. Then copy that substraction into another fits file.
#
# Type: Script
#
# Returns: fits image corresponding to the substraction

__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "2.0.0"
__maintainer__ = "Guadalupe Canas"
__email__ = "gch24@alumnos.unican.es"


import matplotlib.pyplot as plt #Works for making python behave like matlab
import numpy as np #Maths arrays and more
from scipy.special import erfcinv #Inverse Complementary Error Function
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
from astropy.io import fits #Open and Reading Files
import matplotlib.pylab as P #For histograms
import matplotlib

data_image_r = 'lhn1n1_crossmatching_1to2_corrected_leg.fits'
data_image_z = 'lhn1n1_crossmatching_2to1_corrected_leg.fits'

hdulist_data_image_r=fits.open(data_image_r, memmap=False)
hdulist_data_image_z=fits.open(data_image_z, memmap=False)

x_data_image_r=hdulist_data_image_r[0].header['NAXIS1']
y_data_image_r=hdulist_data_image_r[0].header['NAXIS2']

x_data_image_z=hdulist_data_image_z[0].header['NAXIS1']
y_data_image_z=hdulist_data_image_z[0].header['NAXIS2']

data_r =  hdulist_data_image_r[0].data
data_z =  hdulist_data_image_z[0].data

matrix=np.zeros(shape=(y_data_image_r, x_data_image_r))

hdulist_data_image_r.close(closed=True)
hdulist_data_image_z.close(closed=True)

matrix = data_r-data_z

fits.writeto('Substraction.fits', matrix)