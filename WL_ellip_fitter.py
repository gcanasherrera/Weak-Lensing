# Name: WL_ellip_fitter.py
#
# Weak-Lensing "Study of Systematics and Classification of Compact Objects" Program II
#
# Type: python software
#
# Description: From the data obtained thanks to ellipto (which provides a good fitting of stars/galaxies sizes), ellip_fitting.py calculates again the values of ellip_1 and ellip_2, calculates a fitting of those magnitudes f(ellip_) in terms of Legendre Polynomials L(x)L(y) --> f(ellip_)=Pn_m=Cn_m Ln(x)Lm(y) or a normal fitting f(ellip_)=Pn_m=Cn_m x^ny^m, where n=m=4, and plots the results. Then, in makes a iterative fitting rejecting those celestial objects whose residuals with respect to the fitting are bigger than 3*std_deviation(residual).
#
# Returns: FIAT File containing the values of Cn_m in a format suitable to be understood by dlscombine
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "3.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"

# Modifications: introduce the possibility of choosing between a normal Polinomial fitting or Legendre fitting

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma 
import sys
import math
import subprocess
from astropy.modeling import models, fitting #Package for fitting Legeandre Polynomials
import warnings
from mpl_toolkits.mplot3d import Axes3D #Plotting in 3D
import seaborn as sns
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


# 3D Plot function of fitting + experimental data
def ellip_fitting_plt(cont, ellipticity_1, ellipticity_2, fitting_1, fitting_2, x_fit_1, y_fit_1, x_fit_2, y_fit_2):
    sns.set(style="ticks")
    fig = plt.figure(cont)
    ax_1 = fig.add_subplot(211, projection='3d')
    ax_1.plot(x_fit_1, y_fit_1, ellipticity_1, 'k*')
    ax_1.set_xlabel('X ', fontsize=20)
    ax_1.set_ylabel('Y', fontsize=20)
    ax_1.set_zlabel('$e_1$', fontsize=20)
    x_th = np.linspace(0, 12000, 10000)
    y_th = np.linspace(0, 10000, 10000)
    ax_1.plot(x_fit_1, y_fit_1, fitting_1(x_fit_1, y_fit_1), 'b-', label='Legendre fitting')
    ax_2 = fig.add_subplot(212, projection='3d')
    ax_2.plot(x_fit_2, y_fit_2, ellipticity_2, 'k*')
    ax_2.set_xlabel('X', fontsize=20)
    ax_2.set_ylabel('Y', fontsize=20)
    ax_2.set_zlabel('$e_2$', fontsize=20)
    ax_2.plot(x_fit_2, y_fit_2, fitting_2(x_fit_2, y_fit_2), 'b-', label='Legendre fitting')
    return



def ellip_fitting_plt_color(cont, ellipticity_1, ellipticity_2, fitting_1, fitting_2, x_fit_1, y_fit_1, x_fit_2, y_fit_2):


    xv_1, yv_1 = np.meshgrid(x_fit_1,y_fit_1)
    xv_2, yv_2 = np.meshgrid(x_fit_2,y_fit_2)
    levels = MaxNLocator(nbins=15).tick_values(ellipticity_1.min(), ellipticity_1.max())

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    fig, (ax0, ax1) = plt.subplots(nrows=2)

    im = ax0.pcolor(xv_1, yv_1, fitting_1(xv_1,yv_1), cmap='RdBu', norm=norm)
    cbar=fig.colorbar(im, ax=ax0)
    cbar.ax.tick_params(labelsize=30)
    ax0.set_title('$\mathrm{e_1}$ (top) and $\mathrm{e_2}$ (bottom)', fontsize=30)
    ax0.set_xlabel('$x_1/pixels$', fontsize=30, labelpad=10)
    ax0.set_ylabel('$y_1/pixels$', fontsize=30)

    #plt.xticks(color='k', size=30)
    #plt.yticks(color='k', size=30)
    
    
    for tick in ax0.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax0.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)


    im = ax1.pcolormesh(xv_2, yv_2, fitting_2(xv_2, yv_2), cmap='RdBu', norm=norm)
    cbar=fig.colorbar(im, ax=ax1)
    cbar.ax.tick_params(labelsize=30)
    #ax1.set_title('$e_2$', fontsize=20)
    ax1.set_xlabel('$x_2/pixels$', fontsize=30, labelpad=20)
    ax1.set_ylabel('$y_2/pixels$', fontsize=30)

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)


def fit_Legendre(FILE_NAME, catalog):
    
    # Define the Legendre Model in 2D dimensions trough the library ASTROPY
    # The fitting is in the form of Pn_m=Cn_m Ln(x)Lm(y) where the maximum degree is 4, meaning x(degree)*y(degree)<=4.
    # Number of coefficients npar=15
    # Coefficient Matrix:
    # c0_0     c0_1    c0_2     c0_3    c0_4
    # c1_0     c1_1    c1_2     c1_3    0
    # c2_0     c2_1    c2_2     0       0
    # c3_0     c3_1    0        0       0
    # c4_0     0       0        0       0
    #
    legendre_1=models.Legendre2D(x_degree=4, y_degree=4)
    
    # Fixed the coefficients that must be hold 0 during the fitting
    
    legendre_1.c4_1.fixed = True
    legendre_1.c3_2.fixed = True
    legendre_1.c4_2.fixed = True
    legendre_1.c2_3.fixed = True
    legendre_1.c3_3.fixed = True
    legendre_1.c4_3.fixed = True
    legendre_1.c1_4.fixed = True
    legendre_1.c2_4.fixed = True
    legendre_1.c3_4.fixed = True
    legendre_1.c4_4.fixed = True
    
    
    legendre_2=models.Legendre2D(x_degree=4,y_degree=4)
    
    # Fixed the coefficients that must be hold 0 during the fitting
    
    legendre_2.c4_1.fixed = True
    legendre_2.c3_2.fixed = True
    legendre_2.c4_2.fixed = True
    legendre_2.c2_3.fixed = True
    legendre_2.c3_3.fixed = True
    legendre_2.c4_3.fixed = True
    legendre_2.c1_4.fixed = True
    legendre_2.c2_4.fixed = True
    legendre_2.c3_4.fixed = True
    legendre_2.c4_4.fixed = True
    
    #Calculation of ellip_
    
    longitud=len(catalog["theta"])
    ellip_1=np.zeros(longitud)
    ellip_2=np.zeros(longitud)
    ellip_error_1=np.zeros(longitud)
    ellip_error_2=np.zeros(longitud)
    for i in range(0, longitud): #multiply array
        ellip_1[i]=catalog["ellipticity"][i]*math.cos(2*math.radians(catalog["theta"][i]))
        ellip_error_1[i]=ellip_1[i]*(catalog['sigma_e'][i]/catalog["ellipticity"][i])
        #ellip_error_1[i] = catalog['sigma_e'][i]*math.cos(2*math.radians(catalog["theta"][i]))
        ellip_2[i]=catalog["ellipticity"][i]*math.sin(2*math.radians(catalog["theta"][i]))
        ellip_error_2[i]=ellip_2[i]*(catalog['sigma_e'][i]/catalog["ellipticity"][i])
        #ellip_error_2[i] = catalog['sigma_e'][i]*math.sin(2*math.radians(catalog["theta"][i]))

    #Read variables x and y from catalog
    x=catalog["x"]
    y=catalog["y"]

    #Define fitting using ASTROPY library: Simplex algorithm and least squares statistic
    #fit=fitting.SimplexLSQFitter()
    fit=fitting.LevMarLSQFitter()
    #fit=fitting.SLSQPLSQFitter()




    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        #Fit
        fit_1 = fit(legendre_1, x, y, ellip_1, weights=catalog['sigma_e'])
        fit_2 = fit(legendre_2, x, y, ellip_2, weights=catalog['sigma_e'])
    #FIT ELLIP_1
    print 'Values of ellip_1 Legendre fitting: {}'.format(fit_1)
    #FIT ELLIP_2
    print 'Values of ellip_2 Legendre fitting: {}'.format(fit_2)
    
    #Plot first fit. Define contador cont to open new plot windows
    cont=16
    ellip_fitting_plt(cont, ellip_1, ellip_2, fit_1, fit_2, x, y, x, y)
    plt.show(block=False)
    ellip_fitting_plt_color(cont, ellip_1, ellip_2, fit_1, fit_2, x, y, x, y)
    plt.show()
    
    #Second fitting based of boundin residuals values. First define residuals
    residual_1=np.zeros(longitud)
    residual_2=np.zeros(longitud)
    
    #Evaluate the residuals as R=(values_exp-values_theo). Defined as abs value because residual can be either negative or positive but we are just interested in the total magnitude
    for j in range(0, longitud):
        residual_1[j]=np.absolute(ellip_1[j] - fit_1(x[j], y[j]))
        residual_2[j]=np.absolute(ellip_2[j] - fit_2(x[j], y[j]))

    #We use a mask: we hide the values of residual_ that are over a fixed value nsig*stad_dev
    nsig=3
    cut_1=nsig*np.std(residual_1, dtype=np.float64)
    print 'The value of the cut for ellip_1 is equal to cut_1={}'.format(cut_1)
    mask_1=residual_1<cut_1
    cut_2=nsig*np.std(residual_2, dtype=np.float64)
    print 'The value of the cut for ellip_2 is equal to cut_2={}'.format(cut_2)
    mask_2=residual_2<cut_2
    
    #We created masked arrays of x, y and ellip_ according to the mask obtained through mask_
    x_int_1= x[mask_1]
    y_int_1= y[mask_1]
    ellip_1_int= ellip_1[mask_1]
    ellip_error_1_int=catalog['sigma_e'][mask_1]
    x_int_2= x[mask_2]
    y_int_2= y[mask_2]
    ellip_2_int= ellip_2[mask_2]
    ellip_error_2_int=catalog['sigma_e'][mask_2]
    
    #We can re-calculate the fitting
    legendre_1_int=models.Legendre2D(4,4)
    legendre_1_int.c4_1.fixed = True
    legendre_1_int.c3_2.fixed = True
    legendre_1_int.c4_2.fixed = True
    legendre_1_int.c2_3.fixed = True
    legendre_1_int.c3_3.fixed = True
    legendre_1_int.c4_3.fixed = True
    legendre_1_int.c1_4.fixed = True
    legendre_1_int.c2_4.fixed = True
    legendre_1_int.c3_4.fixed = True
    legendre_1_int.c4_4.fixed = True
    legendre_2_int=models.Legendre2D(4,4)
    legendre_2_int.c4_1.fixed = True
    legendre_2_int.c3_2.fixed = True
    legendre_2_int.c4_2.fixed = True
    legendre_2_int.c2_3.fixed = True
    legendre_2_int.c3_3.fixed = True
    legendre_2_int.c4_3.fixed = True
    legendre_2_int.c1_4.fixed = True
    legendre_2_int.c2_4.fixed = True
    legendre_2_int.c3_4.fixed = True
    legendre_2_int.c4_4.fixed = True
    #fit_int=fitting.SimplexLSQFitter()
    fit_int=fitting.LevMarLSQFitter()
    #fit_int=fitting.SLSQPLSQFitter()
    fit_1_int = fit_int(legendre_1_int, x_int_1, y_int_1, ellip_1_int, weights=ellip_error_1_int)
    fit_2_int = fit_int(legendre_2_int, x_int_2, y_int_2, ellip_2_int, weights=ellip_error_2_int)

    
    print 'Values of ellip_1_int: {}'.format(fit_1_int)
    print '\n'
    print '\n'
    #print fit_1_int.fit_info
    print '\n'
    print '\n'
    print 'Values of ellip_2_int: {}'.format(fit_2_int)
    print '\n'
    print '\n'
    
    #Plot iterative fit
    cont=cont+1
    ellip_fitting_plt(cont, ellip_1_int, ellip_2_int, fit_1_int, fit_2_int, x_int_1, y_int_1, x_int_2, y_int_2)
    plt.show()
    ellip_fitting_plt_color(cont, ellip_1_int, ellip_2_int, fit_1_int, fit_2_int, x_int_1, y_int_1, x_int_2, y_int_2)
    plt.show()
    
    #Now we want to write the information into a file that is able to be understood by fiat.c and dlscombine.c.
    
    # read in fit.  The fitfile looks like this (after header):
    # p[0]      p[order+1]   ... ...  p[npar/2-1]
    # p[1]      p[order+2]            0
    # ...       ...                   0
    # ...       p[2*order]            0
    # p[order]  0            0   0    0
    # where xout = p[0]+p[1]*x...+p[order]*x^order+p[order+1]*y
    # +p[order+2]*y*x ... + p[npar/2-1]*y^order
    #
    #
    
    #We made ourselves the FIAT file
    #First: e_1
    fitting_file_ellip= '{}_fit_leg.fcat'.format(FILE_NAME)
    f_1=open(fitting_file_ellip, 'w')
    f_1.write('# fiat 1.0\n')
    f_1.write('# nobj_init = {} / before clipping \n'.format(longitud))
    f_1.write('# nob_fin = {} / after clipping\n'.format(len(x_int_1)))
    f_1.write('# npar = 15 / number of parameters \n')
    f_1.write('# order = 4 / order of polynomial fit \n')
    f_1.write('# nfits = 2 / e_1, e_2 \n')
    f_1.write('# nsig = 3.0 / number of sigma for clipping \n')
    f_1.write('# ttype1 = yorder0 \n')
    f_1.write('# ttype2 = yorder1 \n')
    f_1.write('# ttype3 = yorder2 \n')
    f_1.write('# ttype4 = yorder3 \n')
    f_1.write('# ttype5 = yorder4 \n')
    f_1.write('# e_1 / e_1=e*cos(2*theta)\n')
    f_1.write('# residual_mean = {} \n'. format(np.mean(residual_1, dtype=np.float64)))
    f_1.write('# residual_mean = {} \n'. format(np.std(residual_1, dtype=np.float64)))
    
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c0_0.value, fit_1_int.c0_1.value, fit_1_int.c0_2.value, fit_1_int.c0_3.value, fit_1_int.c0_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c1_0.value, fit_1_int.c1_1.value, fit_1_int.c1_2.value, fit_1_int.c1_3.value, fit_1_int.c1_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c2_0.value, fit_1_int.c2_1.value, fit_1_int.c2_2.value, fit_1_int.c2_3.value, fit_1_int.c2_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c3_0.value, fit_1_int.c3_1.value, fit_1_int.c3_2.value, fit_1_int.c3_3.value, fit_1_int.c3_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c4_0.value, fit_1_int.c4_1.value, fit_1_int.c4_2.value, fit_1_int.c4_3.value, fit_1_int.c4_4.value))
    
    #Second: e_2
    f_1.write('# e_2 / e_2=e*sin(2*theta)\n')
    f_1.write('# residual_mean = {} \n'. format(np.mean(residual_2, dtype=np.float64)))
    f_1.write('# residual_mean = {} \n'. format(np.std(residual_2, dtype=np.float64)))
    
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c0_0.value, fit_2_int.c0_1.value, fit_2_int.c0_2.value, fit_2_int.c0_3.value, fit_2_int.c0_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c1_0.value, fit_2_int.c1_1.value, fit_2_int.c1_2.value, fit_2_int.c1_3.value, fit_2_int.c1_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c2_0.value, fit_2_int.c2_1.value, fit_2_int.c2_2.value, fit_2_int.c2_3.value, fit_2_int.c2_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c3_0.value, fit_2_int.c3_1.value, fit_2_int.c3_2.value, fit_2_int.c3_3.value, fit_2_int.c3_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c4_0.value, fit_2_int.c4_1.value, fit_2_int.c4_2.value, fit_2_int.c4_3.value, fit_2_int.c4_4.value))
    f_1.close()

    return fitting_file_ellip

def fit_Polynomial(FILE_NAME, catalog):
    # Define the Polynomial Model in 2D dimensions trough the library ASTROPY for e_1 and e_2
    # The fitting is in the way of Pn_m=Cn_m x^n y^m where the maximum degree is 4, meaning x(degree)*y(degree)<=4.
    # Number of coefficients npar=15
    # Coefficient Matrix:
    # c0_0     c0_1    c0_2     c0_3    c0_4
    # c1_0     c1_1    c1_2     c1_3    0
    # c2_0     c2_1    c2_2     0       0
    # c3_0     c3_1    0        0       0
    # c4_0     0       0        0       0
    #
    p_init_1 = models.Polynomial2D(degree=4)
    
    p_init_2 = models.Polynomial2D(degree=4)
    
    #Calculation of ellip_
    
    longitud=len(catalog["theta"])
    ellip_1=np.zeros(longitud)
    ellip_2=np.zeros(longitud)
    for i in range(0, longitud): #multiply array
        ellip_1[i]=catalog["ellipticity"][i]*math.cos(2*math.radians(catalog["theta"][i]))
        ellip_2[i]=catalog["ellipticity"][i]*math.sin(2*math.radians(catalog["theta"][i]))

    #Read variables x and y from catalog
    x=catalog["x"]
    y=catalog["y"]

    #Define fitting using ASTROPY library: Simplex algorithm and least squares statistic
    fit=fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        #Fit
        fit_1 = fit(p_init_1, x, y, ellip_1)
        fit_2 = fit(p_init_2, x, y, ellip_2)
    #FIT ELLIP_1
    print 'Values of ellip_1 using Polynomial fitting: {}'.format(fit_1)
    #FIT ELLIP_2
    print 'Values of ellip_2 using Polynomial fitting: {}'.format(fit_2)
    
    #Plot first fit. Define contador cont to open new plot windows
    cont=16
    ellip_fitting_plt(cont, ellip_1, ellip_2, fit_1, fit_2, x, y, x, y)
    plt.show(block=False)
    
    #Second fitting based of boundin residuals values. First define residuals
    residual_1=np.zeros(longitud)
    residual_2=np.zeros(longitud)
    
    #Evaluate the residuals as R=(values_exp-values_theo). Defined as abs value because residual can be either negative or positive but we are just interested in the total magnitude
    for j in range(0, longitud):
        residual_1[j]=np.absolute(ellip_1[j] - fit_1(x[j], y[j]))
        residual_2[j]=np.absolute(ellip_2[j] - fit_2(x[j], y[j]))

    #We use a mask: we hide the values of residual_ that are over a fixed value nsig*stad_dev
    nsig=3
    cut_1=nsig*np.std(residual_1, dtype=np.float64)
    print 'The value of the cut for ellip_1 is equal to cut_1={}'.format(cut_1)
    mask_1=residual_1<cut_1
    cut_2=nsig*np.std(residual_2, dtype=np.float64)
    print 'The value of the cut for ellip_2 is equal to cut_2={}'.format(cut_2)
    mask_2=residual_2<cut_2
    
    #We created masked arrays of x, y and ellip_ according to the mask obtained through mask_
    x_int_1= x[mask_1]
    y_int_1= y[mask_1]
    ellip_1_int= ellip_1[mask_1]
    x_int_2= x[mask_2]
    y_int_2= y[mask_2]
    ellip_2_int= ellip_2[mask_2]
    
    #We can re-calculate the fitting
    p_init_1_int=models.Polynomial2D(degree=4)
    
    p_init_2_int=models.Polynomial2D(degree=4)
    
    fit_int=fitting.LevMarLSQFitter()
    fit_1_int = fit_int(p_init_1_int, x_int_1, y_int_1, ellip_1_int)
    fit_2_int = fit_int(p_init_2_int, x_int_2, y_int_2, ellip_2_int)
    
    print 'Values of ellip_1_int using Polynomial fitting: {}'.format(fit_1_int)
    print 'Values of ellip_2_int using Polynomial fitting: {}'.format(fit_2_int)
    
    #Plot iterative fit
    cont=cont+1
    ellip_fitting_plt(cont, ellip_1_int, ellip_2_int, fit_1_int, fit_2_int, x_int_1, y_int_1, x_int_2, y_int_2)
    plt.show()
    
    #ellip_fitting_plt_color(cont, ellip_1_int, ellip_2_int, fit_1_int, fit_2_int, x_int_1, y_int_1, x_int_2, y_int_2)
    #plt.show()

    
    #Now we want to write the information into a file that is able to be understood by fiat.c and dlscombine.c.
    
    # read in fit.  The fitfile looks like this (after header):
    # p[0]      p[order+1]   ... ...  p[npar/2-1]
    # p[1]      p[order+2]            0
    # ...       ...                   0
    # ...       p[2*order]            0
    # p[order]  0            0   0    0
    # where xout = p[0]+p[1]*x...+p[order]*x^order+p[order+1]*y
    # +p[order+2]*y*x ... + p[npar/2-1]*y^order
    #
    #
    
    #We made ourselves the FIAT file
    #First: e_1
    fitting_file_ellip= '{}_fit_pol.fcat'.format(FILE_NAME)
    f_1=open(fitting_file_ellip, 'w')
    f_1.write('# fiat 1.0\n')
    f_1.write('# nobj_init = {} / before clipping \n'.format(longitud))
    f_1.write('# nob_fin = {} / after clipping\n'.format(len(x_int_1)))
    f_1.write('# npar = 15 / number of parameters \n')
    f_1.write('# order = 4 / order of polynomial fit \n')
    f_1.write('# nfits = 2 / e_1, e_2 \n')
    f_1.write('# nsig = 3.0 / number of sigma for clipping \n')
    f_1.write('# ttype1 = yorder0 \n')
    f_1.write('# ttype2 = yorder1 \n')
    f_1.write('# ttype3 = yorder2 \n')
    f_1.write('# ttype4 = yorder3 \n')
    f_1.write('# ttype5 = yorder4 \n')
    f_1.write('# e_1 / e_1=e*cos(2*theta)\n')
    f_1.write('# residual_mean = {} \n'. format(np.mean(residual_1, dtype=np.float64)))
    f_1.write('# residual_mean = {} \n'. format(np.std(residual_1, dtype=np.float64)))
    
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c0_0.value, fit_1_int.c0_1.value, fit_1_int.c0_2.value, fit_1_int.c0_3.value, fit_1_int.c0_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c1_0.value, fit_1_int.c1_1.value, fit_1_int.c1_2.value, fit_1_int.c1_3.value, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c2_0.value, fit_1_int.c2_1.value, fit_1_int.c2_2.value, 0, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c3_0.value, fit_1_int.c3_1.value, 0, 0, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_1_int.c4_0.value, 0, 0, 0, 0))
    
    #Second: e_2
    f_1.write('# e_2 / e_2=e*sin(2*theta)\n')
    f_1.write('# residual_mean = {} \n'. format(np.mean(residual_2, dtype=np.float64)))
    f_1.write('# residual_mean = {} \n'. format(np.std(residual_2, dtype=np.float64)))
    
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c0_0.value, fit_2_int.c0_1.value, fit_2_int.c0_2.value, fit_2_int.c0_3.value, fit_2_int.c0_4.value))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c1_0.value, fit_2_int.c1_1.value, fit_2_int.c1_2.value, fit_2_int.c1_3.value, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c2_0.value, fit_2_int.c2_1.value, fit_2_int.c2_2.value, 0, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c3_0.value, fit_2_int.c3_1.value, 0, 0, 0))
    f_1.write('%-20s\t%-20s\t%-20s\t%-20s\t%-20s \n'% (fit_2_int.c4_0.value, 0, 0, 0, 0))
    f_1.close()
    return fitting_file_ellip