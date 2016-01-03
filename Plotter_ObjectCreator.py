
from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math

FILTER = 'MH'
PICTURE = 'lhn1n1_2010apr_r_stack_fc_fix'


def m(F, a, b):
    return -2.5*(math.log(F, 10)-math.log(a, 10)-b)



#Read txt file ' w2_53_stack_simulation.txt'


mag_input = np.genfromtxt(('{}_simulation_mag_input_{}.txt').format(PICTURE, FILTER))
mag_output_sex = np.genfromtxt(('{}_simulation_mag_output_sex_{}.txt').format(PICTURE, FILTER))
#mag_output_wayback = np.genfromtxt(('w2_53_stack_simulation_mag_output_wayback_{}.txt').format(FILTER))
mag_output_error_sex = np.genfromtxt(('{}_simulation_mag_output_error_sex_{}.txt').format(PICTURE, FILTER))
#mag_output_error_wayback = np.genfromtxt(('w2_53_stack_simulation_mag_output_error_wayback_{}.txt').format(FILTER))
flux_input = np.genfromtxt(('{}_simulation_flux_input_{}.txt').format(PICTURE, FILTER))
flux_output = np.genfromtxt(('{}_simulation_flux_output_{}.txt').format(PICTURE, FILTER))
flux_output_error = np.genfromtxt(('{}_simulation_flux_output_error_{}.txt').format(PICTURE, FILTER))
flux_output_max = np.genfromtxt(('{}_simulation_flux_output_max_{}.txt').format(PICTURE, FILTER))
flux_output_max_error = np.genfromtxt(('{}_simulation_flux_output_max_error_{}.txt').format(PICTURE, FILTER))
number_lost_objects = np.genfromtxt(('{}_simulation_number_lost_objects_{}.txt').format(PICTURE, FILTER))


param = ["mean_a", "mean_b"]
par = np.genfromtxt(('{}_axis_param_{}.txt').format(PICTURE, FILTER), names=param)

flux_input_iso = []
mag_input_iso = []

#PLOT 1: Number of lost Galaxies vs mag_input
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
#plt.title('Number of lost Galaxies vs mag_input_max')
plt.plot(mag_input, number_lost_objects, 'k-')
plt.xlabel('$m(max)_{input}$', fontsize=18)
plt.ylabel('$n_{lost}$', fontsize=18)
plt.show()


#PLOT 2: flux_out_max vs flux_input_max (Linear Scale)

normalization_th = 2*math.pi*par["mean_a"]*par["mean_b"]

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
#plt.title('flux_out_max vs flux_input (Linear Scale)')
plt.errorbar(flux_input, flux_output_max, flux_output_max_error, 0, fmt='ko')
plt.xlim(0, 1e12)
plt.ylim(0, 1e12)
plt.xlabel('$F(max)_{input}$', fontsize=18)
plt.ylabel('$F(max)_{output}$', fontsize=18)
plt.show()

#PLOT 3: flux_out_iso vs flux_input_iso

ratio_flux=map(truediv, flux_output, flux_input)
normalization_exp = np.mean(ratio_flux)

print '\nTheoretical Normalization: {}'.format(normalization_th)
print '\nExperimental Normalization: {}'.format(ratio_flux)


for i in range(len(flux_input)):
#flux_input_iso.append(ratio_flux[i]*flux_input[i])
    flux_input_iso.append(normalization_th*flux_input[i])

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
#plt.title('flux_out_iso vs flux_input_iso (Linear Scale)')
plt.errorbar(flux_input_iso, flux_output, flux_output_error, 0, fmt='ko')
plt.xlabel('$F(iso)_{input}$', fontsize=18)
plt.ylabel('$F(iso)_{output}$', fontsize=18)
plt.show()


#PLOT 4: mag_output_iso vs mag_input_iso

ratio_mag=map(truediv, mag_output_sex, mag_input)

for i in flux_input_iso:
    mag_input_iso.append(m(i, 100.228418351, 9.99901042564))

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
#plt.title('mag_out_iso vs mag_input_iso (Linear Scale)')
plt.errorbar(mag_input_iso, mag_output_sex, mag_output_error_sex, 0, fmt='ko')
plt.xlabel('$m(iso)_{input}$', fontsize=18)
plt.ylabel('$m(iso)_{output}$', fontsize=18)
plt.show()



print 'END ! ! ! \n'
