
from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math


#Read txt file

mag_input = np.genfromtxt('w2_53_stack_simulation_mag_input.txt')
mag_output_sex = np.genfromtxt('w2_53_stack_simulation_mag_output_sex.txt')
mag_output_wayback = np.genfromtxt('w2_53_stack_simulation_mag_output_wayback.txt')
mag_output_error_sex = np.genfromtxt('w2_53_stack_simulation_mag_output_error_sex.txt')
mag_output_error_wayback = np.genfromtxt('w2_53_stack_simulation_mag_output_error_wayback.txt')
flux_input = np.genfromtxt('w2_53_stack_simulation_flux_input.txt')
flux_output = np.genfromtxt('w2_53_stack_simulation_flux_output.txt')
flux_output_error = np.genfromtxt('w2_53_stack_simulation_flux_output_error.txt')
flux_output_max = np.genfromtxt('w2_53_stack_simulation_flux_output_max.txt')
flux_output_max_error = np.genfromtxt('w2_53_stack_simulation_flux_output_max_error.txt')
number_lost_objects = np.genfromtxt('w2_53_stack_simulation_number_lost_objects.txt')


param = ["mean_a", "mean_b"]
par = np.genfromtxt('w2_53_stack_axis_param.txt', names=param)

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
for i in range (0, len(ratio_flux):
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

for i in range (0, len(ratio_mag):
    mag_input_iso.append(ratio_mag[i]*mag_input[i])

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('mag_out_iso vs mag_input_iso (Linear Scale)')
plt.errorbar(mag_input_iso, mag_output_sex, mag_output_error_sex, 0, fmt='ko')
plt.xlabel('mag_input_iso')
plt.ylabel('mag_output_iso')
plt.show()



print 'END ! ! ! \n'
