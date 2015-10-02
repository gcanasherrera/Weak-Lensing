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

#mag_input = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
mag_input = [6,8,10]
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
    print '\nRound {}\n'.format(mag_input[i])
    print '\nMasking \n'
    simulation.masking_matrix('w2_53_stack.fits')
    simulation.packing_percentage(number_objects = 1000)
    
    
    simulation.out_mag = [0]
    simulation.objectcreator_magnitude(mag_value = mag_input[i], n = 5)
    
    print 'Sextractor'
    sex_caller('w2_53_stack_Simulation_{}.fits'.format(mag_input[i]), 'w2_53_stack_simulation_{}'.format(mag_input[i]))
    catag_simulation = CatalogReader('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    catag_simulation.read('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    simulation.searcher_kdtree(catag.fcat, catag_simulation.fcat, 'w2_stack_53')
    
    print 'The mean value of the output magnitude is {} and the std deviation is {}\n'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag))
    mag_output.append(np.mean(simulation.out_mag))
    mag_output_error.append(np.std(simulation.out_mag))
    print

np.savetxt('mag_ouput2.txt', mag_output, delimiter=',')

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, mag_output, 0, mag_output_error, fmt='ko', label='Different Magnitudes')
plt.xlabel('Mag_input')
plt.ylabel('Mag_output')
plt.show()




print 'END ! ! ! \n'
