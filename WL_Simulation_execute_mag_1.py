# Name: WL_ObjectCreator_execute_mag.py
#
# Bachelor Disertation Program VI
#
# Type: Script
#
# Description: simulation in magnitude --> creates diverse celestial objects with the same magnitude mag_input (remind this is the maximum value) and evaluate the discrepancies with Source Extractor
#
# Returns: text file with values obtaining from the simulation
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain GNU"
__version__ = "2.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"

from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math



#Define arrays for plotting
mag_input = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

mag_output_sex = []
mag_output_wayback = []
mag_output_error_sex = []
mag_output_error_wayback = []
flux_input = []
flux_output = []
flux_output_error = []
flux_output_max = []
flux_output_max_error = []
number_lost_objects = []

FILTER = 'Mexican_hat'


#Read catalog
catag = CatalogReader('w2_53_stack.fcat')
catag.read()

#Create object for simulation
simulation = ObjectCreator(catag.fcat)

print '\nMasking . . . {}\n'

simulation.masking_matrix('w2_53_stack.fits')

#For-loop for 1 to 30 mag
for mag in mag_input:
    
    print '\nRound {}\n'.format(mag)
    simulation.packing_percentage(number_objects = 4000)
    simulation.out_mag = []
    simulation.posible_obj_distances = []
    simulation.posible_obj_index = []
    simulation.out_flux = []
    simulation.out_mag_after_transf = []
    simulation.out_flux_max = []
    
    
    
    print '\nSimulation\n'
    
    simulation.objectcreator_magnitude(mag_value = mag, n = 5)
    
    print '\nSextractor\n'
    sex_caller('w2_53_stack_simulation_{}.fits'.format(mag), 'w2_53_stack_simulation_{}'.format(mag))
    catag_simulation = CatalogReader('w2_53_stack_simulation_{}.cat'.format(mag))
    catag_simulation.read()
    
    
    print '\nSearcher\n'
    simulation.searcher_kdtree(catag.fcat, catag_simulation.fcat, 'w2_stack_53')
    
    print 'The mean value of the output sextractor magnitude is {}\nThe std deviation is {}\nThe mean value of the output wayback magnitude is {}\nThe std deviation (wayback magnitude) is {}\nThe mean value of the output flux is {}\nThe std deviation (flux) is {}\n'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag), np.mean(simulation.out_mag_after_transf), np.std(simulation.out_mag_after_transf), np.mean(simulation.out_flux), np.std(simulation.out_flux))
    
    #Save values into arrays
    
    mag_output_sex.append(np.mean(simulation.out_mag))
    mag_output_wayback.append(np.mean(simulation.out_mag_after_transf))
    mag_output_error_sex.append(np.std(simulation.out_mag))
    mag_output_error_wayback.append(np.std(simulation.out_mag_after_transf))
    flux_output.append(np.mean(simulation.out_flux))
    flux_output_error.append(np.std(simulation.out_flux))
    flux_output_max.append(np.mean(simulation.out_flux_max))
    flux_output_max_error.append(np.std(simulation.out_flux_max))
    np.savetxt('simulation_index_{}.txt'.format(mag), simulation.posible_obj_index)


#Getting flux_max_input from mag_input

fitting_flux = MagnitudeExponential()

for mag in mag_input:
    flux_input.append(fitting_flux.f(mag, a=simulation.parameter_a,b=simulation.parameter_b))

#SAVE DATA in txt file

np.savetxt('w2_53_stack_simulation_mag_input_{}.txt'.format(FILTER), mag_input)
np.savetxt('w2_53_stack_simulation_mag_output_sex_{}.txt'.format(FILTER), mag_output_sex)
np.savetxt('w2_53_stack_simulation_mag_output_wayback_{}.txt'.format(FILTER), mag_output_wayback)
np.savetxt('w2_53_stack_simulation_mag_output_error_sex_{}.txt'.format(FILTER), mag_output_error_sex)
np.savetxt('w2_53_stack_simulation_mag_output_error_wayback_{}.txt'.format(FILTER), mag_output_error_wayback)
np.savetxt('w2_53_stack_simulation_flux_input_{}.txt'.format(FILTER), flux_input)
np.savetxt('w2_53_stack_simulation_flux_output_{}.txt'.format(FILTER), flux_output)
np.savetxt('w2_53_stack_simulation_flux_output_error_{}.txt'.format(FILTER), flux_output_error)
np.savetxt('w2_53_stack_simulation_flux_output_max_{}.txt'.format(FILTER), flux_output_max)
np.savetxt('w2_53_stack_simulation_flux_output_max_error_{}.txt'.format(FILTER), flux_output_max_error)
np.savetxt('w2_53_stack_simulation_number_lost_objects_{}.txt'.format(FILTER), simulation.lost_objects)


np.savetxt('w2_53_stack_axis_param_{}.txt'.format(FILTER), (simulation.mean_a, simulation.mean_b))

