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
from operator import truediv
from Class_ObjectCreator import magnitude_exponential
import math


catag = CatalogReader('w2_53_stack.fcat')
catag.read('w2_53_stack.fcat')

simulation = ObjectCreator(catag.fcat)
#simulation.general_histograms(simulation.fcat)

mag_input = [5]
#mag_input = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
#mag_input = [19, 20, 21, 22, 23, 24, 25, 26, 27]
#mag_input = [6]
mag_output_sex = []
mag_output_wayback = []
mag_output_error_sex = []
mag_output_error_wayback = []
flux_input = []
flux_output = []
flux_output_error = []
flux_output_max = []
flux_output_max_error = []



for i in range (0, len(mag_input)):
    print '\nRound {}\n'.format(mag_input[i])
    print '\nMasking \n'
    simulation.masking_matrix('w2_53_stack.fits')
    simulation.packing_percentage(number_objects = 10)
    
    simulation.out_mag = []
    simulation.posible_obj_distances = []
    simulation.posible_obj_index = []
    simulation.out_flux = []
    simulation.out_mag_after_transf = []
    
    print '\nSimulation\n'
    
    simulation.objectcreator_magnitude(mag_value = mag_input[i], n = 5)
    
    print '\nSextractor\n'
    sex_caller('w2_53_stack_Simulation_{}.fits'.format(mag_input[i]), 'w2_53_stack_simulation_{}'.format(mag_input[i]))
    catag_simulation = CatalogReader('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    catag_simulation.read('w2_53_stack_simulation_{}.cat'.format(mag_input[i]))
    
    
    print '\nSearcher\n'
    simulation.searcher_kdtree_try(catag.fcat, catag_simulation.fcat, 'w2_stack_53')
    
    print 'The mean value of the output sextractor magnitude is {}\nThe std deviation is {}\nThe mean value of the output wayback magnitude is {}\nThe std deviation (wayback magnitude) is {}\nThe mean value of the output flux is {}\nThe std deviation (flux) is {}\n'.format(np.mean(simulation.out_mag), np.std(simulation.out_mag), np.mean(simulation.out_mag_after_transf), np.std(simulation.out_mag_after_transf), np.mean(simulation.out_flux), np.std(simulation.out_flux))
    
    mag_output_sex.append(np.mean(simulation.out_mag))
    mag_output_wayback.append(np.mean(simulation.out_mag_after_transf))
    mag_output_error_sex.append(np.std(simulation.out_mag))
    mag_output_error_wayback.append(np.std(simulation.out_mag_after_transf))
    flux_output.append(np.mean(simulation.out_flux))
    flux_output_error.append(np.std(simulation.out_flux))
    flux_output_max.append(np.mean(simulation.out_flux_max))
    flux_output_max_error.append(np.std(simulation.out_flux_max))
    np.savetxt('simulation_index_{}.txt'.format(i), simulation.posible_obj_index)


#PLOT
np.savetxt('simulation_position.txt', simulation.posible_obj_index, delimiter=',')
np.savetxt('mag_ouput_sex.txt', mag_output_sex, delimiter=',')
np.savetxt('mag_ouput_wayback.txt', mag_output_wayback, delimiter=',')
np.savetxt('mag_ouput_sex_error.txt', mag_output_error_sex, delimiter=',')
np.savetxt('mag_ouput_wayback_error.txt', mag_output_error_wayback, delimiter=',')
np.savetxt('flux_output.txt', flux_output, delimiter=',')
np.savetxt('flux_output_error.txt', flux_output_error, delimiter=',')


#PLOT 1: Mag_output_sex vs mag_input
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, mag_output_sex, mag_output_error_sex, 0, fmt='ko', label='Source Extractor')
plt.xlabel('Mag_input')
plt.ylabel('Mag_output_sex')
plt.show(block=False)

#PLOT 2: Mag_output_wayback vs mag_input

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, mag_output_wayback, mag_output_error_wayback, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('Mag_input')
plt.ylabel('Mag_output_wayback')
plt.show(block=False)


#PLOT 3: Number of lost Galaxies vs mag_input

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.plot(mag_input, simulation.lost_objects)
plt.xlabel('Mag_input')
plt.ylabel('Number Lost Objects')
plt.show(block=False)

#PLOT 4: mag_input/mag_out_sex vs mag_input
ratio_mag_sex=map(truediv, mag_input, mag_output_sex)
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, ratio_mag_sex, mag_output_error_sex, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('Mag_input')
plt.ylabel('ratio mag_input/mag_output_sex')
plt.show(block=False)

#PLOT 5: mag_input/mag_out_wayback vs mag_input
ratio_mag_wayback=map(truediv, mag_input, mag_output_wayback)
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(mag_input, ratio_mag_sex, mag_output_error_sex, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('Mag_input')
plt.ylabel('ratio mag_input/mag_output_waveback')
plt.show(block=False)


#PLOT 6: flux_input/flux_out vs flux_input

fitting_flux = magnitude_exponential()

for mag in mag_input:
    flux_input.append(fitting_flux.f(mag, a=112.188243402,b=9.95005045126))

np.savetxt('flux_input.txt', flux_input, delimiter=',')

ratio_flux=map(truediv, flux_input, flux_output)
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(flux_input, ratio_flux, flux_output_error, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('flux_input')
plt.ylabel('ratio flux_input/flux_output')
plt.show(block=False)


#PLOT 7: flux_out_max vs flux_input
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(flux_input, flux_output_max, flux_output_max_error, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('flux_input_max')
plt.ylabel('flux_output_max')
plt.show(block=False)

#PLOT 8: flux_out_max vs flux_input
normalization_th = 2*math.pi*simulation.mean_b*simulation.mean_a
normalization_exp = np.mean(map(truediv, flux_output, flux_input))

print '\nTheoretical Normalization: {}'.format(normalization_th)
print '\nExperimental Normalization: {}'.format(normalization_exp)

flux_input_iso_corrected =  [x * normalization_exp for x in flux_input]


sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.errorbar(flux_input_iso_corrected, flux_output, flux_output_error, 0, fmt='ko', label='Wayback Relation')
plt.xlabel('flux_input_iso')
plt.ylabel('flux_output_iso')
plt.show(block=False)



print 'END ! ! ! \n'