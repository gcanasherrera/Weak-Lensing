# Name: Class_CrossMatching.py
#
# Bachelor Disertation Program X
#
# Type: python class
#
# Content: 2 Classes, 1 constructor,
#
# Description: General Class destinated to get a cross matching between two catalogs from different filtered pics


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
from astropy.io import fits #Open and Reading FITS Files usign astropy
import random #pseudo-random generator
import seaborn as sns #Improvements for statistical-plots
from WL_Utils import sex_caller
from Class_CrossMatching import CrossMatching


#Get first catalogs
sex_caller('lhn1n1_2010apr_r_stack_fc_fix.fits'.format(mag), 'lhn1n1_2010apr_r_stack_fc_fix')
sex_caller('lhn1n1_2010dec_z_stack_fc_fix.fits'.format(mag), 'lhn1n1_2010dec_z_stack_fc_fix')

#Read catalogs
catag_r = CatalogReader('lhn1n1_2010apr_r_stack_fc_fix.fcat')
catag_r.read()

catag_z = CatalogReader('lhn1n1_2010dec_z_stack_fc_fix.fcat')
catag_z.read()

#Create object for cross-matching

crossmatching = CrossMatching(catag_r.fcat, catag_z.fcat)
crossmatching.kdtree(n=2)
crossmatching.catalog_writter('lhn1n1_crossmatching_1to2.fcat', compare = '1to2'):



