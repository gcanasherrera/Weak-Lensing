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
class CrossMatching(object):
    
    
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


    def positions_array_creator(self):
        self.positions_cat_1 = zip(self.x_1.ravel(), self.y_1.ravel())
        self.positions_cat_2 = zip(self.x_2.ravel(), self.y_2.ravel())

    def kdtree(self, n=1):
        self.positions_array_creator()
        tree_1 = spatial.KDTree(self.positions_cat_1)
        self.posible_obj_distances_1, self.posible_obj_index_1 = tree_1.query(self.positions_cat_2, distance_upper_bound = n)

        tree_2 = spatial.KDTree(self.positions_cat_2)
        self.posible_obj_distances_2, self.posible_obj_index_2 = tree_2.query(self.positions_cat_1, distance_upper_bound = n)

    def catalog_writter(self, FILE_NAME, compare = ''):
        catalog_def= '{}.fcat'.format(FILE_NAME)
        f_1=open(catalog_def, 'w')
        f_1.write('# fiat 1.0\n')
        f_1.write('# written by Class_CrossMatching.py\n')
        f_1.write('# TTYPE1 = NUMBER / Running object number \n')
        f_1.write('# TTYPE2 = FLUX_ISO / Isophotal flux \n')
        f_1.write('# TTYPE3 = FLUXERR_ISO / RMS error for isophotal flux   \n')
        f_1.write('# TTYPE4 = MAG_ISO / Isophotal magnitude    \n')
        f_1.write('# TTYPE5 = MAGERR_ISO / RMS error for isophotal magnitude   \n')
        f_1.write('# TTYPE6 = MAG_APER_1 / Fixed aperture magnitude       \n')
        f_1.write('# TTYPE7 = MAGERR_APER_1 / RMS error  for fixed aperture mag.   \n')
        f_1.write('# TTYPE8 = MAG / Best of MAG_AUTO and MAG_ISOCOR   \n')
        f_1.write('# TTYPE9 = MAGERR / RMS error for MAG_BEST     \n')
        f_1.write('# TTYPE10 = FLUX_MAX / Peak flux above background    \n')
        f_1.write('# TTYPE11 = ISOAREA / Isophotal area above Analysis threshold    \n')
        f_1.write('# TTYPE12 = X / Object position along x   \n')
        f_1.write('# TTYPE13 = Y / Object position along y    \n')
        f_1.write('# TTYPE14 = ra / Right ascension of barycenter (J2000)  \n')
        f_1.write('# TTYPE15 = dec / Declination of barycenter (J2000)   \n')
        f_1.write('# TTYPE16 = ixx / Variance along x     \n')
        f_1.write('# TTYPE17 = iyy / Variance along y    \n')
        f_1.write('# TTYPE18 = ixy / Covariance between x and y   \n')
        f_1.write('# TTYPE19 = ixxWIN / Windowed variance along x    \n')
        f_1.write('# TTYPE20 = iyyWIN / Windowed variance along y    \n')
        f_1.write('# TTYPE21 = ixyWIN / Windowed covariance between x and y    \n')
        f_1.write('# TTYPE22 = A / Profile RMS along major axis    \n')
        f_1.write('# TTYPE23 = B / Profile RMS along minor axis    \n')
        f_1.write('# TTYPE24 = THETA / Position angle (CCW/x)    \n')
        f_1.write('# TTYPE25 = ELONGATION / A_IMAGE/B_IMAGE    \n')
        f_1.write('# TTYPE26 = ELLIPTICITY / 1 - B_IMAGE/A_IMAGE    \n')
        f_1.write('# TTYPE27 = FWHM / FWHM assuming a gaussian core    \n')
        f_1.write('# TTYPE28 = FLAGS / Extraction flags    \n')
        f_1.write('# TTYPE29 = CLASS_STAR / S/G classifier output     \n')

        cont = 1
        number_of_lost_1 = 0
        number_of_lost_2 = 0
        
        #define array of catalog names
        
        names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]

        if compare == '1to2':

            for index in self.posible_obj_index_1:
                if index < len(self.positions_cat_1):
                
                    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (cont, self.catalog_1[names[1]][index], self.catalog_1[names[2]][index], self.catalog_1[names[3]][index], self.catalog_1[names[4]][index], self.catalog_1[names[5]][index], self.catalog_1[names[6]][index], self.catalog_1[names[7]][index], self.catalog_1[names[8]][index], self.catalog_1[names[9]][index], self.catalog_1[names[10]][index], self.catalog_1[names[11]][index], self.catalog_1[names[12]][index], self.catalog_1[names[13]][index], self.catalog_1[names[14]][index], self.catalog_1[names[15]][index], self.catalog_1[names[16]][index], self.catalog_1[names[17]][index], self.catalog_1[names[18]][index], self.catalog_1[names[19]][index], self.catalog_1[names[20]][index], self.catalog_1[names[21]][index], self.catalog_1[names[22]][index], self.catalog_1[names[23]][index], self.catalog_1[names[24]][index], self.catalog_1[names[25]][index], self.catalog_1[names[26]][index], self.catalog_1[names[27]][index], self.catalog_1[names[28]][index]))
                
                    cont = cont + 1

                elif index == len(self.positions_cat_1):
                    number_of_lost_1 = number_of_lost_1 + 1
    
            print 'The number of lost objects from matching of 1 with respect of 2 is {}\n'.format(number_of_lost_1)
            print 'Length of the Picture 1 original catalog {}\n'.format(len(self.catalog_1[names[0]]))
            print 'Length of the Cross-Matching Catalog {}\n'.format(cont)


#if compare == '2to1':
#
#            for index in self.posible_obj_index_2:
#                if index < len(self.positions_cat_2):
                    
                    #                    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (cont, self.catalog_2[names[1]], self.catalog_2[names[2]], self.catalog_2[names[3]], self.catalog_2[names[4]], self.catalog_2[names[5]], self.catalog_2[names[6]], self.catalog_2[names[7]], self.catalog_2[names[8]], self.catalog_2[names[9]], self.catalog_2[names[10]], self.catalog_2[names[11]], self.catalog_2[names[12]], self.catalog_2[names[13]], self.catalog_2[names[14]], self.catalog_2[names[15]], self.catalog_2[names[16]], self.catalog_2[names[17]], self.catalog_2[names[18]], self.catalog_2[names[19]], self.catalog_2[names[20]], self.catalog_2[names[21]], self.catalog_2[names[22]], self.catalog_2[names[23]], self.catalog_2[names[24]], self.catalog_2[names[25]], self.catalog_2[names[26]], self.catalog_2[names[27]], self.catalog_2[names[28]]))
                        
            #  cont = cont + 1
                            
                            #    elif index == len(self.positions_cat_2):
                            #number_of_lost_2 = number_of_lost_2 + 1

# print 'The number of lost objects from matching of 2 with respect of 1 is {}\n'.format(number_of_lost_2)
#    print 'Length of the Picture 2 original catalog {}\n'.format(len(self.catalog_2[names[0]]))
#      print 'Length of the Cross-Matching Catalog {}\n'.format(cont)

