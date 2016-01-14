import seaborn as sns 
import numpy as np
import math
from scipy.stats import norm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt


def gaussian2D(x, y, x_0, y_0, a, b, M):
    return M * math.exp(-0.5*((x-x_0)*(x-x_0)/(a*a) + (y-y_0)*(y-y_0)/(b*b)))


mu_x, sigma_x = 1532.27929344, 1.3 # mean and standard deviation
s_x = np.random.normal(mu_x, sigma_x, 1000000)

mu_y, sigma_y = 1532.27929344, 0.936 # mean and standard deviation
s_y = np.random.normal(mu_y, sigma_y, 1000000)

#sns.jointplot(s_x, s_y, kind="kde")

sns.set(style="ticks")
plt.hist2d(s_x, s_y, bins=100, norm=LogNorm())
plt.colorbar()
plt.show()
sns.jointplot(s_x, s_y, kind="kde", size=7, space=0)