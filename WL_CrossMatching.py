# Name: WL_CrossMatching.py
#
# Bachelor Disertation Program II
#
# Type: python script
#
# Content: two objects from the class CatalogReader, and 1 object from the class CrossMatching
#
# Description: Main Script to perform the Cross Matching algorithm


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
from Class_CatalogReader import CatalogReader


#Get first catalogs
#sex_caller('lhn1n1_2010apr_r_stack_fc_fix.fits', 'lhn1n1_2010apr_r_stack_fc_fix')
#sex_caller('lhn1n1_2010dec_z_stack_fc_fix.fits', 'lhn1n1_2010dec_z_stack_fc_fix')

#Read catalogs
catag_r = CatalogReader('lhn1n1_2010apr_r_stack_fc_fix.cat')
catag_r.read()

catag_z = CatalogReader('lhn1n1_2010dec_z_stack_fc_fix.cat')
catag_z.read()

#Give value to the cross-matching radius

r=3

#Create object for cross-matching
crossmatching = CrossMatching(catag_r.fcat, catag_z.fcat)
crossmatching.kdtree(n=r*1e-06)
crossmatching.catalog_writter('lhn1n1_crossmatching_1to2', compare = '1to2')
print '\n'
crossmatching.catalog_writter('lhn1n1_crossmatching_2to1', compare = '2to1')

if crossmatching.cont1to2<crossmatching.cont2to1:
    catag_final_1 = CatalogReader('lhn1n1_crossmatching_1to2.fcat')
    catag_final_1.read()
    catag_final_2 = CatalogReader('lhn1n1_crossmatching_2to1.fcat')
    catag_final_2.read()
    crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
    crossmatching_final.kdtree(n=r*1e-06)
    crossmatching_final.catalog_writter('lhn1n1_crossmatching_final', compare = '2to1')

if crossmatching.cont1to2>crossmatching.cont2to1:
    catag_final_1 = CatalogReader('lhn1n1_crossmatching_1to2.fcat')
    catag_final_1.read()
    catag_final_2 = CatalogReader('lhn1n1_crossmatching_2to1.fcat')
    catag_final_2.read()
    crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
    crossmatching_final.kdtree(n=r*1e-06)
    crossmatching_final.catalog_writter('lhn1n1_crossmatching_final', compare = '1to2')

if crossmatching.cont1to2==crossmatching.cont2to1:
    catag_final_1 = CatalogReader('lhn1n1_crossmatching_1to2.fcat')
    catag_final_1.read()
    catag_final_2 = CatalogReader('lhn1n1_crossmatching_2to1.fcat')
    catag_final_2.read()
    crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
    crossmatching_final.kdtree(n=r*1e-06)
    crossmatching_final.catalog_writter('lhn1n1_crossmatching_final', compare = '1to2')