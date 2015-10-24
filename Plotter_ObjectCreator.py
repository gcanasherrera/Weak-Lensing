
from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math

FILTER = Gaussian

#Read txt file ' w2_53_stack_simulation.txt'


mag_input = np.genfromtxt(('w2_53_stack_simulation_mag_input_{}.txt').format(FILTER))
mag_output_sex = np.genfromtxt(('w2_53_stack_simulation_mag_output_sex_{}.txt').format(FILTER))
mag_output_wayback = np.genfromtxt(('w2_53_stack_simulation_mag_output_wayback_{}.txt').format(FILTER))
mag_output_error_sex = np.genfromtxt(('w2_53_stack_simulation_mag_output_error_sex_{}.txt').format(FILTER))
mag_output_error_wayback = np.genfromtxt(('w2_53_stack_simulation_mag_output_error_wayback_{}.txt').format(FILTER))
flux_input = np.genfromtxt(('w2_53_stack_simulation_flux_input_{}.txt').format(FILTER))
flux_output = np.genfromtxt(('w2_53_stack_simulation_flux_output_{}.txt').format(FILTER))
flux_output_error = np.genfromtxt(('w2_53_stack_simulation_flux_output_error_{}.txt').format(FILTER))
flux_output_max = np.genfromtxt(('w2_53_stack_simulation_flux_output_max_{}.txt').format(FILTER))
flux_output_max_error = np.genfromtxt(('w2_53_stack_simulation_flux_output_max_error_{}.txt').format(FILTER))
number_lost_objects = np.genfromtxt(('w2_53_stack_simulation_number_lost_objects_{}.txt').format(FILTER))


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

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('flux_out_max vs flux_input (Linear Scale)')
plt.errorbar(flux_input, flux_output_max, flux_output_max_error, 0, fmt='ko')
plt.xlabel('flux_input_max')
plt.ylabel('flux_output_max')
plt.show()

#PLOT 3: flux_out_iso vs flux_input_iso

ratio_flux=map(truediv, flux_output, flux_input)
normalization_exp = np.mean(ratio_flux)

print '\nTheoretical Normalization: {}'.format(normalization_th)
print '\nExperimental Normalization: {}'.format(ratio_flux)


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
