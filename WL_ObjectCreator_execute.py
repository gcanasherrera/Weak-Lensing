# Name: WL_ObjectCreator_execute1.py
#
# Bachelor Disertation Program VI
#
# Type: Script
#
# Description: Gets a
#
# Returns: FITS image including the objects we created
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "1.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"

from Class_ObjectCreator import ObjectCreator
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller

catag = CatalogReader()
catag.read('w2_53_stack.fcat')

simulation = ObjectCreator(catag.fcat)
#simulation.general_histograms(simulation.fcat)
simulation.masking_matrix('w2_53_stack.fits')
simulation.packing_percentage(eta=0.05)
simulation.objectcreator_mean(mag_value = 15, n = 5)


sex_caller('w2_53_stack_Simulation_15.fits', 'w2_53_stack_simulation')
catag_simulation = CatalogReader()
catag_simulation.read('w2_53_stack_simulation.cat')

simulation.searcher(catag.fcat, catag_simulation.fcat)
print 'The mean value of the output magnitude is {} and the std deviation is {}'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag))


print 'END ! ! ! \n'
