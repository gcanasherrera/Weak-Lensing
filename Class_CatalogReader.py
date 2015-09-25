# Name: ObjectCreator.py
#
# Bachelor Disertation Program VI
#
# Type: Script
#
# Description:
#
# Returns: FITS image including the objects we created
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "1.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"


import numpy as np #Maths arrays and more, matlab-type vectors/arrays
import subprocess #calling to the terminal
import sys #Strings inputs

class CatalogReader(object):
    
    def __init__(self):
        
        self.fcat =[]
        
    def transform(self, catalog):
        catalog_fiat = 'fiat_{}'.format(catalog)
        transform_into_fiat='perl sex2fiat.pl {}>{}'.format(catalog, catalog_fiat)
        subprocess.call(transform_into_fiat, shell=True)
    
    def read(self, catalog):
        
        suffix = '.fcat'
        
        names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
        
        if catalog.endswith(suffix) == False:
            transform(catalog)
            self.fcat = np.genfromtxt(catalog_fiat, names=names)

        else:
            self.fcat = np.genfromtxt(catalog, names=names)