import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys
import math
import subprocess
from astropy.modeling import models, fitting #Package for fitting Legeandre Polynomials
import warnings
from mpl_toolkits.mplot3d import Axes3D #Plotting in 3D
#import seaborn as sns

from Class_CatalogReader import CatalogReader

stars_1 = CatalogReader('lhn1n1_2010apr_r_stack_fc_fix_corrected_stars.fcat')
stars_2 = CatalogReader('lhn1n1_2010apr_r_stack_fc_fix_stars.fcat')

stars_1.read()
stars_2.read()

#sns.set_style("white")
#sns.set_style("ticks")


plt.figure()
for line in stars_1.fcat:
    if line['ellipticity']<0.005:
        plt.plot(line['x'], line['y'], 'bo')
    elif line['ellipticity']>0.005:
        plt.plot(line['x'], line['y'], 'k|', markersize=20)



plt.xlabel('x', fontsize=20)
plt.ylabel('y', fontsize=20)
plt.show(block=False)



plt.figure()
for line in stars_2.fcat:
    if line['ellipticity']<0.005:
        plt.plot(line['x'], line['y'], 'bo')
    elif line['ellipticity']>0.005:
        plt.plot(line['x'], line['y'], 'k|', markersize=20)


plt.xlabel('x', fontsize=20)
plt.ylabel('y', fontsize=20)
plt.show()

