
from Class_ObjectCreator import ObjectCreator, MagnitudeExponential
from Class_CatalogReader import CatalogReader
from WL_Utils import sex_caller
import numpy as np #Maths arrays and more
import matplotlib.pyplot as plt #Plot Libraries
import seaborn as sns #Improvements for statistical-plots
from operator import truediv
import math


#Read txt file

variables = ["mag_input", "mag_ouput_sex", "mag_ouput_sex_error", "mag_ouput_wayback", "mag_ouput_wayback_error", "flux_output", "flux_output_error", "flux_output_max", "flux_output_max_error", "flux_input", "number_lost_objects"]

param = ["mean_a, mean_b"]

file = np.genfromtxt('w2_53_stack_simulation_mag.txt', names=variables)
par = np.genfromtxt(''w2_53_stack_axis_param.txt', names=param)


#PLOT 1: Number of lost Galaxies vs mag_input
sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('Number of lost Galaxies vs mag_input_max')
plt.plot(file["mag_input"], file["number_lost_objects"], 'ko')
plt.xlabel('Mag_input_max')
plt.ylabel('Number Lost Objects')
plt.show()


#PLOT 2: flux_out_max vs flux_input_max (Linear Scale)

normalization_th = 2*math.pi*par["mean_a"]*par["mean_b"]
normalization_exp = np.mean(map(truediv, file["flux_output], file[flux_input"]))

print '\nTheoretical Normalization: {}'.format(normalization_th)
print '\nExperimental Normalization: {}'.format(normalization_exp)

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('flux_out_max vs flux_input (Linear Scale)')
plt.errorbar(file["flux_input"], file["flux_output_max"], file["flux_output_max_error"], 0, fmt='ko')
plt.xlabel('flux_input_max')
plt.ylabel('flux_output_max')
plt.show()

#PLOT 3: flux_out_iso vs flux_input_iso

ratio_flux=map(truediv, file["flux_output"], file["flux_input"])
flux_input_iso =  [x * ratio_flux for x in file["flux_input"]]

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('flux_out_iso vs flux_input_iso (Linear Scale)')
plt.errorbar(file["flux_input_iso"], file["flux_output"], file["flux_output_error"], 0, fmt='ko')
plt.xlabel('flux_input_iso')
plt.ylabel('flux_output_iso')
plt.show()


#PLOT 4: mag_output_iso vs mag_input_iso

ratio_mag=map(truediv, mag_output_sex, mag_input)
mag_input_iso =  [x * ratio_mag for x in mag_input]

sns.set(style="white", palette="muted", color_codes=True)
plt.figure()
plt.title('mag_out_iso vs mag_input_iso (Linear Scale)')
plt.errorbar(file["mag_input_iso"], file["mag_output_sex"], file["mag_output_error_sex"], 0, fmt='ko')
plt.xlabel('mag_input_iso')
plt.ylabel('mag_output_iso')
plt.show()



print 'END ! ! ! \n'