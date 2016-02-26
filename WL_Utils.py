# Name: WL_utils.py
#
# Bachelor Disertation Program A.II
#
# Description: Python File that contains all basic functions required by WL_script.py.
#
# For further information about the functions here defined the users are kindly referred to the DOCUMENTATION file attached in this package and also to the coments of each method
#


__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "4.0.0"
__maintainer__ = "Guadalupe Canas Herrera"
__email__ = "gch24@alumnos.unican.es"


import matplotlib.pyplot as plt #Works for making python behave like matlab
#import sextutils as sex #Program used to read the original catalog
import numpy as np #Maths arrays and more
import numpy.ma as ma #Masking arrays
import sys #Strings inputs
import math #mathematical functions
import subprocess #calling to the terminal
from astropy.modeling import models, fitting #Package for fitting Legendre Polynomials
import warnings #Advices
from mpl_toolkits.mplot3d import Axes3D #Plotting in 3D
import WL_ellip_fitter as ellip_fit #Ellipticity fitting
import seaborn as sns


type_cat = ".cat"
type_fcat = ".fcat"
type_good = "_good.fcat"
type_galaxies = "_galaxies.fcat"
type_stars = "_stars.fcat"
type_ellipto_galaxies = "_ellipto_galaxies.fcat"
type_ellipto_stars = "_ellipto_stars.fcat"
type_shapes_galaxies = "_shapes_galaxies.fcat"
type_shapes_stars = "_shapes_stars.fcat"


#
# Function: works to call Source Extractor to obtain celestial objects. You pass the name of the image and the name of the file (catalogue name without ending)
#
def sex_caller(fits, FILE_NAME):
    catalog_name = '{}.cat'.format(FILE_NAME)
    sex= 'sex -CATALOG_NAME {} -CHECKIMAGE_NAME {}_trial.fits {}'.format(catalog_name, FILE_NAME, fits)
    subprocess.call(sex, shell=True)
    return catalog_name

#
# Function: works to call Source Extractor to obtain celestial objects. You pass the name of the image_corrected and the name of the file (catalogue name without ending)
#
def sex_caller_corrected(fits, FILE_NAME):
    catalog_name = '{}_corrected.cat'.format(FILE_NAME)
    sex= 'sex -CATALOG_NAME {} -CHECKIMAGE_NAME {}_trial_corrected.fits {}'.format(catalog_name, FILE_NAME, fits)
    subprocess.call(sex, shell=True)
    return catalog_name


#
# Function: works to call ellipto to re-obtain the sizes of celestial objects. You pass the in - catalog name and the out catalog name
#
def ellipto_caller(catalog_object, fits, catalog_ellipto):
    ellipto= './ellipto -m 512 {} {}>{}'.format(catalog_object, fits, catalog_ellipto)
    subprocess.call(ellipto, shell=True)
    return

#
# Function: works to dlscombine_pol to correct PSF anisotropies. You pass the specfile file and the output name of the corrected photo
#
def dlscombine_pol_caller(specfile, fits_corrected): #Call ellipto when you pass the catalog name
    dlscombine_pol= './dlscombine_pol {} {}'.format(specfile, corrected_photo)
    subprocess.call(dlscombine_pol, shell=True)
    return corrected_photo

#
# Function: works to dlscombine_leg to correct PSF anisotropies. You pass the specfile file and the output name of the corrected photo
#
def dlscombine_leg_caller(specfile, fits_corrected): #Call ellipto when you pass the catalog name
    dlscombine_leg= './dlscombine_leg {} {}'.format(specfile, corrected_photo)
    subprocess.call(dlscombine_leg, shell=True)
    return corrected_photo

#
# Function: works to call ds9 to open FITS files. You pass the name of the image
#
def ds9_caller(photo):
    ds9= './ds9 {}'.format(photo)
    subprocess.call(ds9, shell=True)
    return

#
# Function: Plots whatever st. array is passed as attributes (whatever you want)
#
def plotter(catalog, magnitude_1, magnitude_2, int, legend_1, legend_2):
    sns.set(style="ticks")
    plt.figure(int)
    plt.plot(catalog[magnitude_1], catalog[magnitude_2], 'k.')
    plt.xlabel(legend_1, labelpad=20, fontsize=20)
    plt.ylabel(legend_2, fontsize=20)
    return

#
# Function: Calculate the ellipticity vector components present in the catalog that is passed as attribute. Then, it plots.
#
def ellipticity(catalog, int):
    sns.set(style="ticks")
    leng=len(catalog["theta"])
    ellip_minus=np.zeros(leng)
    ellip_plus=np.zeros(leng)
    for i in range(0, leng): #multiply array
        ellip_minus[i]=catalog["ellipticity"][i]*math.cos(2*math.radians(catalog["theta"][i]))
        ellip_plus[i]=catalog["ellipticity"][i]*math.sin(2*math.radians(catalog["theta"][i]))
    plt.figure(int)
    plt.plot(ellip_minus, ellip_plus, 'k*')
    plt.xlabel('$e_1$', labelpad=20, fontsize=44)
    plt.ylabel('$e_2$', fontsize=44)
    #plt.title('$ellipticity_{minus}$ ' 'Vs.' '$ellipticity_{plus}$' 'for $mag_{iso}$={} and $mag_{iso}$={}'.format(np.amin(catalog['mag_iso']), np.amax(catalog['mag_iso'])))
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.xticks(color='k', size=30)
    plt.yticks(color='k', size=30)
    return

#
# Function: creates a file able to be understood by DLSCombine for the application of the Kernel
#
def specfile(image, ellipfit, FILE_NAME):
    if ellipfit.endswith('_leg.fcat'):
        specfile='{}_specfile_leg.tmp'.format(FILE_NAME)
        f_1=open(specfile, 'w')
        f_1.write('# DLSMAKE = {}\n'.format('CatalogPlotter makes it'))
        f_1.write('# CTYPE1 = {}\n'.format('RA---TAN'))
        f_1.write('# CTYPE2 = {}\n'.format('DEC--TAN'))
        f_1.write('# CRVAL1 = 163.3239204316\n')
        f_1.write('# CRVAL2 = 57.7986803985 \n')
        f_1.write('# CRPIX1 = 5633. \n')
        f_1.write('# CRPIX2 = 4385.\n')
        f_1.write('# CD1_1 = -5.61111E-05\n')
        f_1.write('# CD1_2 = 0.\n')
        f_1.write('# CD2_1 = 0.\n')
        f_1.write('# CD2_2 = 5.61111E-05\n')
        #1.0 unidad 0.00 algo de xy y la ultima deberia ser 60 que es como el peso que se le da a cada cosa pero que para nosotros un numero positivo es suficiente.
        #Ponemos un guion porque queremos que aplique la unidad
        f_1.write('{}\t -\t {}\t 1.0\t 0.000\t 60.000\t'.format(image, ellipfit))
        f_1.close()
    elif ellipfit.endswith('_pol.fcat'):
        specfile='{}_specfile_pol.tmp'.format(FILE_NAME)
        f_1=open(specfile, 'w')
        f_1.write('# DLSMAKE = {}\n'.format('CatalogPlotter makes it'))
        f_1.write('# CTYPE1 = {}\n'.format('RA---TAN'))
        f_1.write('# CTYPE2 = {}\n'.format('DEC--TAN'))
        f_1.write('# CRVAL1 = 163.3239204316\n')
        f_1.write('# CRVAL2 = 57.7986803985 \n')
        f_1.write('# CRPIX1 = 5633. \n')
        f_1.write('# CRPIX2 = 4385.\n')
        f_1.write('# CD1_1 = -5.61111E-05\n')
        f_1.write('# CD1_2 = 0.\n')
        f_1.write('# CD2_1 = 0.\n')
        f_1.write('# CD2_2 = 5.61111E-05\n')
        #1.0 unidad 0.00 algo de xy y la ultima deberia ser 60 que es como el peso que se le da a cada cosa pero que para nosotros un numero positivo es suficiente.
        #Ponemos un guion porque queremos que aplique la unidad
        f_1.write('{}\t -\t {}\t 1.0\t 0.000\t 60.000\t'.format(image, ellipfit))
        f_1.close()
    return specfile


#
# Function: creates a file able to be understood by DLSCombine for the application of the Kernel. Defined for the central pixel (RA and DEC) of the R-band
#

def specfile_r(image, ellipfit, FILE_NAME):
    specfile='{}_specfile_leg.tmp'.format(FILE_NAME)
    f_1=open(specfile, 'w')
    f_1.write('# DLSMAKE = {}\n'.format('CatalogPlotter makes it'))
    f_1.write('# CTYPE1 = {}\n'.format('RA---TAN'))
    f_1.write('# CTYPE2 = {}\n'.format('DEC--TAN'))
    f_1.write('# CRVAL1 = 163.3328929836\n')
    f_1.write('# CRVAL2 = 57.7969314494\n')
    f_1.write('# CRPIX1 = 5571. \n')
    f_1.write('# CRPIX2 = 4388.\n')
    f_1.write('# CD1_1 = -5.61111E-05\n')
    f_1.write('# CD1_2 = 0.\n')
    f_1.write('# CD2_1 = 0.\n')
    f_1.write('# CD2_2 = 5.61111E-05\n')
    #1.0 unidad 0.00 algo de xy y la ultima deberia ser 60 que es como el peso que se le da a cada cosa pero que para nosotros un numero positivo es suficiente.
    #Ponemos un guion porque queremos que aplique la unidad
    f_1.write('{}\t -\t {}\t 1.0\t 0.000\t 60.000\t'.format(image, ellipfit))
    f_1.close()
    return specfile

#
# Function: creates a file able to be understood by DLSCombine for the application of the Kernel. Defined for the central pixel (RA and DEC) of the Z-band
#

def specfile_z(image, ellipfit, FILE_NAME):
    specfile='{}_specfile_leg.tmp'.format(FILE_NAME)
    f_1=open(specfile, 'w')
    f_1.write('# DLSMAKE = {}\n'.format('CatalogPlotter makes it'))
    f_1.write('# CTYPE1 = {}\n'.format('RA---TAN'))
    f_1.write('# CTYPE2 = {}\n'.format('DEC--TAN'))
    f_1.write('# CRVAL1 = 163.3239204316\n')
    f_1.write('# CRVAL2 = 57.7986803985 \n')
    f_1.write('# CRPIX1 = 5633. \n')
    f_1.write('# CRPIX2 = 4385.\n')
    f_1.write('# CD1_1 = -5.61111E-05\n')
    f_1.write('# CD1_2 = 0.\n')
    f_1.write('# CD2_1 = 0.\n')
    f_1.write('# CD2_2 = 5.61111E-05\n')
    #1.0 unidad 0.00 algo de xy y la ultima deberia ser 60 que es como el peso que se le da a cada cosa pero que para nosotros un numero positivo es suficiente.
    #Ponemos un guion porque queremos que aplique la unidad
    f_1.write('{}\t -\t {}\t 1.0\t 0.000\t 60.000\t'.format(image, ellipfit))
    f_1.close()
    return specfile



#
# Function: filter only stars (two criteria at the same time)
#
def stars_maker(catalog, file_name):
    print("Let's obtain only a FIAT catalog that contains stars. We need to bound. Have a look to the FWHM vs Mag_ISO plot")
    mag_iso_min_stars=float(raw_input("Enter the minimum value for mag_iso: "))
    mag_iso_max_stars=float(raw_input("Enter the maximum value for mag_iso: "))
    FWHM_min_stars=float(raw_input("Enter the minimum value for FWHM: "))
    FWHM_max_stars=float(raw_input("Enter the maximum value for FWHM: "))
    catalog_name_stars= '{}{}'.format(file_name, type_stars)
    #Creamos un string para que lo ponga en la terminal
    terminal_stars= 'perl fiatfilter.pl "MAG_ISO>{} && MAG_ISO<{} && FWHM>{} && FWHM<{} && CLASS_STAR>0.9 && FLUX_ISO<3000000" {}>{}'.format(mag_iso_min_stars, mag_iso_max_stars, FWHM_min_stars, FWHM_max_stars, catalog, catalog_name_stars)
    subprocess.call(terminal_stars, shell=True)
    return catalog_name_stars, FWHM_max_stars

#
# Function: filter only galaxies performing a linear fit (two criteria at the same time)
#
def galaxies_maker(catalog, file_name, FWHM_max_stars):
    print("Let's obtain only a FIAT catalog that contains galaxies. We need to bound. Have a look to the FWHM vs Mag_ISO plot")
    print("")
    print("First, I'm going to perform a linear fit. Tell me the values of mag_iso")
    mag_iso_min_galaxies=float(raw_input("Enter the minimum value for mag_iso: "))
    mag_iso_max_galaxies=float(raw_input("Enter the maximum value for mag_iso: "))
    catalog_name_fit='{}_fit{}'.format(file_name, type_galaxies)
    #Creamos un string para que lo ponga en la terminal
    terminal_fit= 'perl fiatfilter.pl -v "MAG_ISO>{} && MAG_ISO<{}" {}>{}'.format(mag_iso_min_galaxies, mag_iso_max_galaxies, catalog, catalog_name_fit)
    subprocess.call(terminal_fit, shell=True)
    names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
    fcat_fit = np.genfromtxt(catalog_name_fit, names=names)
    fit=np.polyfit(fcat_fit['mag_iso'], fcat_fit['FWHM'], 1)
    #Save in variables the values of the fitting
    m=fit[0]
    n=fit[1]
    print 'The value of the y-intercep n={} and the value of the slope m={}'.format(n,m)
    # Once you have the values of the fitting we can obtain the catalog of galaxies
    catalog_name_galaxies= '{}{}'.format(file_name, type_galaxies)
    terminal_galaxies= 'perl fiatfilter.pl -v "FWHM>{}*MAG_ISO+{} && FWHM>{} && CLASS_STAR<0.1 && FLUX_ISO<3000000" {}>{}'.format(m, n, FWHM_max_stars, catalog, catalog_name_galaxies)
    subprocess.call(terminal_galaxies, shell=True)
    return catalog_name_galaxies
