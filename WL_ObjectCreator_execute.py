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
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots


catag = CatalogReader('w2_53_stack.fcat')
catag.read('w2_53_stack.fcat')

simulation = ObjectCreator(catag.fcat)
#simulation.general_histograms(simulation.fcat)
simulation.masking_matrix('w2_53_stack.fits')
simulation.packing_percentage(eta=0.05)

mag_input = [7, 8, 9, 10]
mag_output = []

mag_output_error = []

#simulation.objectcreator_magnitude(mag_value = 2, n = 5)
#sex_caller('w2_53_stack_Simulation_{}.fits'.format(mag_input[i]), 'w2_53_stack_simulation_{}'.format(mag_input[i]))
#catag_simulation = CatalogReader('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
#catag_simulation.read('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
#simulation.searcher_dic(catag.fcat, catag_simulation.fcat)
#print 'The mean value of the output magnitude is {} and the std deviation is {}'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag))
#mag_output.append(simulation.out_mag)

for i in range (0, len(mag_input)):
    print mag_input[i]
    simulation.objectcreator_magnitude(mag_value = mag_input[i], n = 5)
    sex_caller('w2_53_stack_Simulation_{}.fits'.format(mag_input[i]), 'w2_53_stack_simulation_{}'.format(mag_input[i]))
    catag_simulation = CatalogReader('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    catag_simulation.read('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    simulation.searcher_dic(catag.fcat, catag_simulation.fcat)
    print 'The mean value of the output magnitude is {} and the std deviation is {}'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag))
    mag_output.append(np.mean(simulation.out_mag))
    mag_output_error.append(np.std(simulation.out_mag))


sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, mag_output, np.std(mag_input), np.std(mag_output), fmt='k.', label='Different Magnitudes')
plt.xlabel('Mag_input')
plt.ylabel('Mag_output')
plt.show()




print 'END ! ! ! \n'
