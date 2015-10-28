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
from pymodelfit import FunctionModel1DAuto #Create own model of fitting
from scipy import spatial #KDTREE altorithm



"""
    General Class that contents several simulation methods destinated to get a cross-matching between two catalogs
    
"""
class ObjectCreator(object):
    
    
    """
        Constructor that defines attribute of future class-objects
        
    """
    
    def __init__(self, catalog_1, catalog_2):  #define attributes of future class-objects
        # catalog_1 and catalog_2 are objects from Class_CatalogReader.py
        self.catalog_1 = catalog_1
        self.catalog_2 = catalog_2
        self.x_1 = catalog_1['x']
        self.x_2 = catalog_2['x']
        self.y_1 = catalog_1['y']
        self.y_2 = catalog_2['y']
        self.positions_cat_1 = []
        self.positions_cat_2 = []
        self.posible_obj_distances_1 = []
        self.posible_obj_index_1 = []
        self.posible_obj_distances_2 = []
        self.posible_obj_index_2 = []


    def positions_array_creator():
        self.positions_cat_1 = zip(self.x_1.ravel(), self.y_1.ravel())
        self.positions_cat_2 = zip(self.x_2.ravel(), self.y_2.ravel())

    def kdtree(n):
        self.positions_array_creator()
        tree_1 = spatial.KDTree(self.positions_cat_1)
        self.posible_obj_distances_1, self.posible_obj_index_1 = tree_1.query(self.positions_cat_2, distance_upper_bound = n)

        tree_2 = spatial.KDTree(self.positions_cat_2)
        self.posible_obj_distances_2, self.posible_obj_index_2 = tree_2.query(self.positions_cat_1, distance_upper_bound = n)

    def catalog_writter():
        fitting_file_ellip= '{}_fit_pol.fcat'.format(FILE_NAME)
        f_1=open(fitting_file_ellip, 'w')
        f_1.write('# fiat 1.0\n')
        f_1.write('# nobj_init = {} / before clipping \n'.format(longitud))
        f_1.write('# nob_fin = {} / after clipping\n'.format(len(x_int_1)))
        f_1.write('# npar = 15 / number of parameters \n')
        f_1.write('# order = 4 / order of polynomial fit \n')
        f_1.write('# nfits = 2 / e_1, e_2 \n')
        f_1.write('# nsig = 3.0 / number of sigma for clipping \n')
        f_1.write('# ttype1 = yorder0 \n')
        f_1.write('# ttype2 = yorder1 \n')
        f_1.write('# ttype3 = yorder2 \n')
        f_1.write('# ttype4 = yorder3 \n')
        f_1.write('# ttype5 = yorder4 \n')
        f_1.write('# e_1 / e_1=e*cos(2*theta)\n')
        f_1.write('# residual_mean = {} \n'. format(np.mean(residual_1, dtype=np.float64)))
        f_1.write('# residual_mean = {} \n'. format(np.std(residual_1, dtype=np.float64)))
    
        f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c0_0.value, fit_1_int.c0_1.value, fit_1_int.c0_2.value, fit_1_int.c0_3.value, fit_1_int.c0_4.value))





















