# Name: ObjectCreator.py
#
# Bachelor Disertation Program A.I
#
# Type: Python Class
#
# Content: 1 class, 1 constructor, 2 methods
#
# Description: Read Source Extractor and FIAT catalogues transforming information into structured arrays from numpy
#
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "2.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"


import numpy as np #Maths arrays and more, matlab-type vectors/arrays
import subprocess #calling to the terminal
import sys #Strings inputs

"""
    General Class that contents several methods destinated to read catalogues
    
"""

class CatalogReader(object):
    
    """
        Constructor that reads the name of the catalogue and split the name into the name itself and the ending
        
    """
    
    def __init__(self, catalog):
        
        self.catalog = catalog
        self.BEFORE_NAME = self.catalog.find('.')
        self.FILE_NAME = self.catalog[:self.BEFORE_NAME]
        
        self.fcat =[]
        self.catalog_fiat=''
    
    """
        Method that transform the Source Extractor catalogue into FIAT catalogue
        
    """
        
    def transform(self):
        self.catalog_fiat = 'fiat_{}.fcat'.format(self.FILE_NAME)
        transform_into_fiat='perl sex2fiat.pl {}>{}'.format(self.catalog, self.catalog_fiat)
        subprocess.call(transform_into_fiat, shell=True)
    
    
    """
        Method that reads the catalogue and transforms it into numpy arrays
        
    """
    
    def read(self):
        
        suffix = '.fcat'
        names = ["x", "y", "mag_iso", "median", "ixx", "iyy", "ixy", "a_input", "b_input", "theta", "ellipticity", "errcode", "sigsky", "size", "flux", "mean_rho_4th", "sigma_e", "wander"]
        #names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
        
        #For older version of Source Extractor
        
        names_cat = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
        
        if self.catalog.endswith(suffix) == False:
            self.transform()
            self.fcat = np.genfromtxt(self.catalog_fiat, names=names)

        else:
            self.fcat = np.genfromtxt(self.catalog, names=names)


