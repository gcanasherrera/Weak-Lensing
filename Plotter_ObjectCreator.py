
from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math


#Read txt file ' w2_53_stack_simulation.txt'


mag_input_1 = np.genfromtxt('w2_53_stack_simulation_mag_input_1.txt')
mag_output_sex_1 = np.genfromtxt('w2_53_stack_simulation_mag_output_sex_1.txt')
mag_output_wayback_1 = np.genfromtxt('w2_53_stack_simulation_mag_output_wayback_1.txt')
mag_output_error_sex_1 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_sex_1.txt')
mag_output_error_wayback_1 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_wayback_1.txt')
flux_input_1 = np.genfromtxt('w2_53_stack_simulation_flux_input_1.txt')
flux_output_1 = np.genfromtxt('w2_53_stack_simulation_flux_output_1.txt')
flux_output_error_1 = np.genfromtxt('w2_53_stack_simulation_flux_output_error_1.txt')
flux_output_max_1 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_1.txt')
flux_output_max_error_1 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_error_1.txt')
number_lost_objects_1 = np.genfromtxt('w2_53_stack_simulation_number_lost_objects_1.txt')



mag_input_2 = np.genfromtxt('w2_53_stack_simulation_mag_input_2.txt')
mag_output_sex_2 = np.genfromtxt('w2_53_stack_simulation_mag_output_sex_2.txt')
mag_output_wayback_2 = np.genfromtxt('w2_53_stack_simulation_mag_output_wayback_2.txt')
mag_output_error_sex_2 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_sex_2.txt')
mag_output_error_wayback_2 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_wayback_2.txt')
flux_input_2 = np.genfromtxt('w2_53_stack_simulation_flux_input_2.txt')
flux_output_2 = np.genfromtxt('w2_53_stack_simulation_flux_output_2.txt')
flux_output_error_2 = np.genfromtxt('w2_53_stack_simulation_flux_output_error_2.txt')
flux_output_max_2 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_2.txt')
flux_output_max_error_2 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_error_2.txt')
number_lost_objects_2 = np.genfromtxt('w2_53_stack_simulation_number_lost_objects_2.txt')


mag_input_3 = np.genfromtxt('w2_53_stack_simulation_mag_input_3.txt')
mag_output_sex_3 = np.genfromtxt('w2_53_stack_simulation_mag_output_sex_3.txt')
mag_output_wayback_3 = np.genfromtxt('w2_53_stack_simulation_mag_output_wayback_3.txt')
mag_output_error_sex_3 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_sex_3.txt')
mag_output_error_wayback_3 = np.genfromtxt('w2_53_stack_simulation_mag_output_error_wayback_3.txt')
flux_input_3 = np.genfromtxt('w2_53_stack_simulation_flux_input_3.txt')
flux_output_3 = np.genfromtxt('w2_53_stack_simulation_flux_output_3.txt')
flux_output_error_3 = np.genfromtxt('w2_53_stack_simulation_flux_output_error_3.txt')
flux_output_max_3 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_3.txt')
flux_output_max_error_3 = np.genfromtxt('w2_53_stack_simulation_flux_output_max_error_3.txt')
number_lost_objects_3 = np.genfromtxt('w2_53_stack_simulation_number_lost_objects_3.txt')



mag_input = np.append(mag_input_1, mag_input_2, mag_input_3)
mag_output_sex = np.append(mag_output_sex_1, mag_output_sex_2, mag_output_sex_3)
mag_output_wayback = np.append(mag_output_wayback_1, mag_output_wayback_2, mag_output_wayback_3)
mag_output_error_sex = np.append(mag_output_error_sex_1, mag_output_error_sex_2, mag_output_error_sex_3)
mag_output_error_wayback = np.append(mag_output_error_wayback_1, mag_output_error_wayback_2, mag_output_error_wayback_3)
flux_input = np.append(flux_input_1, flux_input_2, flux_input_3)
flux_output = np.append(flux_output_1, flux_output_2, flux_output_3)
flux_output_error = np.append(flux_output_error_1, flux_output_error_2, flux_output_error_3)
flux_output_max = np.append(flux_output_max_1, flux_output_max_2, flux_output_max_3)
flux_output_max_error = np.append(flux_output_max_error_1, flux_output_max_error_2, flux_output_max_error_3)
number_lost_objects = np.append(number_lost_objects_1, number_lost_objects_2, number_lost_objects_3)




param = ["mean_a", "mean_b"]
par = np.genfromtxt('w2_53_stack_axis_param_1.txt', names=param)

flux_input_iso = []
mag_input_iso = []

#PLOT 1: Number of lost Galaxies vs mag_input
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('Number of lost Galaxies vs mag_input_max')
plt.plot(mag_input, number_lost_objects, 'ko')
plt.xlabel('Mag_input_max')
plt.ylabel('Number Lost Objects')
plt.show()


#PLOT 2: flux_out_max vs flux_input_max (Linear Scale)

normalization_th = 2*math.pi*par["mean_a"]*par["mean_b"]
normalization_exp = np.mean(map(truediv, flux_output, flux_input))

print '\nTheoretical Normalization: {}'.format(normalization_th)
print '\nExperimental Normalization: {}'.format(normalization_exp)

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('flux_out_max vs flux_input (Linear Scale)')
plt.errorbar(flux_input, flux_output_max, flux_output_max_error, 0, fmt='ko')
plt.xlabel('flux_input_max')
plt.ylabel('flux_output_max')
plt.show()

#PLOT 3: flux_out_iso vs flux_input_iso

ratio_flux=map(truediv, flux_output, flux_input)
for i in range (0, len(ratio_flux)):
    flux_input_iso.append(ratio_flux[i]*flux_input[i])

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('flux_out_iso vs flux_input_iso (Linear Scale)')
plt.errorbar(flux_input_iso, flux_output, flux_output_error, 0, fmt='ko')
plt.xlabel('flux_input_iso')
plt.ylabel('flux_output_iso')
plt.show()


#PLOT 4: mag_output_iso vs mag_input_iso

ratio_mag=map(truediv, mag_output_sex, mag_input)

for i in range (0, len(ratio_mag)):
    mag_input_iso.append(ratio_mag[i]*mag_input[i])

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('mag_out_iso vs mag_input_iso (Linear Scale)')
plt.errorbar(mag_input_iso, mag_output_sex, mag_output_error_sex, 0, fmt='ko')
plt.xlabel('mag_input_iso')
plt.ylabel('mag_output_iso')
plt.show()



print 'END ! ! ! \n'
