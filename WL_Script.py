# Name: WL_Script.py
#
# Weak-Lensing "Study of Systematics and Classification of Compact Objects" Program I
#
# Type: python script
#
# Description: Central script that develops the whole process of reading images, filtering into galaxies and stars, correcting sizes and shapes, correcting PSF annisotropies, and re-classify compact objects into galaxies to obtain a final catalogue
#
# Returns: FITS image - mass-density map
#          Catalogs
#          Plots
#          FITS image - trial from Source Extractor
#

__author__ = "Guadalupe Canas Herrera"
__copyright__ = "Copyright (C) 2015 G. Canas Herrera"
__license__ = "Public Domain"
__version__ = "4.0.0"
__maintainer__ = "Guadalupe Canas"
__email__ = "gch24@alumnos.unican.es"


# Improvements: more automatic ---> only needs the name of the picture, the catalogue (in case you have it) and the BAND you want to analize
# Old CatalogPlotter3.py has been splitted in two: WL_Script.py and WL_Utils.py
# Also call: WL_utils.py, WL_filter_mag_gal.py - WL_ellip_fitter.py (written by Guadalupe Canas Herrera)
# Also call 2: Source Extractor (by Emmanuel Bertin V2.3.2), sex2fiat (by DAVID WITTMAN v1.2), fiatfilter (by DAVID WITTMAN v1.2), ellipto (by DAVID WITTMAN v1.2), dlscombine (by DAVID WITTMAN v1.2 and modified by GUADALUPE CANAS)

#
# DLSCOMBINE CORRECTS PSF: it has a dependence in fiat.c, fiat.h, dlscombine_utils.c, dlscombine.c, dlscombine.h
# Guadalupe Canas Herrera modified fiat.c, dlscombine_utils.c, dlscombine.h
#


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
from WL_Utils import sex_caller, sex_caller_corrected, ellipto_caller, dlscombine_pol_caller, dlscombine_leg_caller, ds9_caller, plotter, ellipticity, specfile, stars_maker, galaxies_maker, specfile_r, specfile_z
from WL_filter_mag_gal import filter_mag #Filtering final catalog of galaxies a function of magnitudes and call fiatmap
import seaborn as sns
import matplotlib.pylab as P #histograms
from Class_CrossMatching import CrossMatching
from Class_CatalogReader import CatalogReader


###############################  BEGIN SCRIPT   ###############################

# (1): We define the ending of the input/output files

type_fits = ".fits"
type_cat = ".cat"
type_fcat = ".fcat"
type_good = "_good.fcat"
type_galaxies = "_galaxies.fcat"
type_stars = "_stars.fcat"
type_ellipto_galaxies = "_ellipto_galaxies.fcat"
type_ellipto_stars = "_ellipto_stars.fcat"
type_shapes_galaxies = "_shapes_galaxies.fcat"
type_shapes_stars = "_shapes_stars.fcat"
type_match = "_match.fcat"



def main():
    
    sns.set(style="white", palette="muted", color_codes=True)
    
    print("Welcome to the Weak-Lensing Script, here to help you analizing Subaru images  in search of galaxy clusters")
    print("")
    
    array_file_name = []
    
    # (1): Ask the number of image that did the cross-matching process.
    
    question = int(raw_input("Please, tell me how many pictures did the cross-matching: "))
    cont = 0
    BEFORE_NAME = ''
    FILE_NAME = ''
    #print FILE_NAME
    FILE_NAME_CORRECTED= ''
    
    while cont < question:
        
        # (2): We need to read the image and band. We ask in screen the image of the region of the sky.
    
        filter =raw_input("Introduce the name of the filter: ")
        fits = raw_input("Please, introduce the name of the fits image you want to read or directly the catalogue: ")
        

        #Save the name of the .fits and .cat in a string:
        
        BEFORE_NAME = fits.find('.')
        FILE_NAME = fits[:BEFORE_NAME]
        #print FILE_NAME
        FILE_NAME_CORRECTED='{}_corrected'.format(FILE_NAME)
        
        
        if fits.endswith(type_fits):
        
            #(3) STEP: Call Source Extractor
            
            print("Let me call Source Extractor (called sex by friends). It will obtain the celestial objects. When it finishes I will show you the trial image")
            print("")
            
            catalog_name = sex_caller(fits, FILE_NAME)
            #Show results of trial.fits
            #subprocess.call('./ds9 {}_trial.fits'.format(FILE_NAME), shell=True)
            
            
            #(4): Transform Source Extractor catalog into FIAT FORMAT
            print("I'm transforming the catalog into a FIAT 1.0 format")
            print("")
            
            catalog_name_fiat= '{}.fcat'.format(FILE_NAME)
            transform_into_fiat='perl sex2fiat.pl {}>{}'.format(catalog_name, catalog_name_fiat)
            subprocess.call(transform_into_fiat, shell=True)
        
        if fits.endswith(type_fcat):
            catalog_name_fiat = fits
            fits = raw_input("Please, introduce the name of the fits image: ")
        
        #(5): Read the FIAT Catalog
        FWHM_max_stars=0
        names = ["number", "flux_iso", "fluxerr_iso", "mag_iso", "magger_iso", "mag_aper_1", "magerr_aper_1", "mag", "magger", "flux_max", "isoarea", "x", "y", "ra", "dec", "ixx", "iyy", "ixy", "ixxWIN", "iyyWIN", "ixyWIN", "A", "B", "theta", "enlogation", "ellipticity", "FWHM", "flags", "class_star"]
        fcat = np.genfromtxt(catalog_name_fiat, names=names)
        P.figure()
        P.hist(fcat['class_star'], 50, normed=1, histtype='stepfilled')
        P.show()
        
        #Let's fix the ellipcity + and - for all celestial objects

        
        #(6): plot FWHM vs mag_iso
        print("I'm ploting MAG_ISO vs. FWHM")
        magnitude1='mag_iso'
        magnitude2='FWHM'
        plotter(fcat, magnitude1, magnitude2, 2, '$mag(iso)$', '$FWHM/pixels$')
        plt.show()
        print("Do you want to fix axis limits? Please answer with y or n")
        answer=raw_input()
        if answer== "y":
            xmin=float(raw_input("X min: "))
            xmax=float(raw_input("X max: "))
            ymin=float(raw_input("Y min: "))
            ymax=float(raw_input("Y max: "))
            #Fix limits
            plotter(catalog_name, magnitude1, magnitude2, 3)
            plt.xlim(xmin,xmax)
            plt.ylim(ymin,ymax)
            plt.show(block=False)
        elif answer == "n":
            plt.show(block=False)
        else:
            plt.show(block=False)

        # (7): Obtaining a GOOD CATALOG without blank spaces and filter saturate objects
        print("This catalog is not the good one. I'll show you why")
        print("")
        magnitude_x="x"
        magnitude_y="y"
        plotter(fcat, magnitude_x, magnitude_y, 4, '$x/pixels$', '$y/pixels$')
        plt.show(block=False)
        print("Please, introduce the values you prefer to bound x and y")
        xmin_good=float(raw_input("X min: "))
        xmax_good=float(raw_input("X max: "))
        ymin_good=float(raw_input("Y min: "))
        ymax_good=float(raw_input("Y max: "))
        catalog_name_good= '{}{}'.format(FILE_NAME, type_good)
        terminal_good= 'perl fiatfilter.pl "x>{} && x<{} && y>{} && y<{} && FLUX_ISO<3000000" {}>{}'.format(xmin_good, xmax_good, ymin_good, ymax_good, catalog_name_fiat, catalog_name_good)
        subprocess.call(terminal_good, shell=True)
        print("Wait a moment, I'm showing you the results in a sec")
        fcat_good = np.genfromtxt(catalog_name_good, names=names)
        print np.amax(fcat_good['flux_iso'])
        plotter(fcat_good, 'x', 'y', 5, '$x/pixels$', '$y/pixels$')
        plt.show(block=False)
        ellipticity(fcat_good, 1)
        plt.show(block=False)
        
        plotter(fcat_good, magnitude1, magnitude2, 2, '$mag(iso)$', '$FWHM/pixels$')
        plt.show(block=False)


        #(8.1.): Creating STARS CATALOG
        print("Let's obtain only a FIAT catalog that contains stars. We need to bound. Have a look to the FWHM vs Mag_ISO plot")
        mag_iso_min_stars=float(raw_input("Enter the minimum value for mag_iso: "))
        mag_iso_max_stars=float(raw_input("Enter the maximum value for mag_iso: "))
        FWHM_min_stars=float(raw_input("Enter the minimum value for FWHM: "))
        FWHM_max_stars=float(raw_input("Enter the maximum value for FWHM: "))
        catalog_name_stars= '{}{}'.format(FILE_NAME, type_stars)
        #Creamos un string para que lo ponga en la terminal
        terminal_stars= 'perl fiatfilter.pl "MAG_ISO>{} && MAG_ISO<{} && FWHM>{} && FWHM<{} && CLASS_STAR>0.9 && FLUX_ISO<3000000" {}>{}'.format(mag_iso_min_stars, mag_iso_max_stars, FWHM_min_stars, FWHM_max_stars, catalog_name_good, catalog_name_stars)
        subprocess.call(terminal_stars, shell=True)
        fcat_stars=np.genfromtxt(catalog_name_stars, names=names)
        ellipticity(fcat_stars, 6)
        plt.show(block=False)

        #(8.2.): Checking STARS CATALOG with Source Extractor Neural Network Output
        P.figure()
        P.hist(fcat_stars['class_star'], 50, normed=1, histtype='stepfilled')
        P.show(block=False)



        #(9.1.): Creating GALAXIES CATALOG
        print("Let's obtain only a FIAT catalog that contains galaxies. We need to bound. Have a look to the FWHM vs Mag_ISO plot")
        print("")
        print("First, I'm going to perform a linear fit. Tell me the values of mag_iso")
        mag_iso_min_galaxies=float(raw_input("Enter the minimum value for mag_iso: "))
        mag_iso_max_galaxies=float(raw_input("Enter the maximum value for mag_iso: "))
        catalog_name_fit='{}_fit{}'.format(FILE_NAME, type_galaxies)
        #Creamos un string para que lo ponga en la terminal
        terminal_fit= 'perl fiatfilter.pl -v "MAG_ISO>{} && MAG_ISO<{}" {}>{}'.format(mag_iso_min_galaxies, mag_iso_max_galaxies, catalog_name_good, catalog_name_fit)
        subprocess.call(terminal_fit, shell=True)
        fcat_fit = np.genfromtxt(catalog_name_fit, names=names)
        fit=np.polyfit(fcat_fit['mag_iso'], fcat_fit['FWHM'], 1)
        #Save in variables the values of the fitting
        m=fit[0]
        n=fit[1]
        print 'The value of the y-intercep n={} and the value of the slope m={}'.format(n,m)
        # Once you have the values of the fitting we can obtain the catalog of galaxies
        catalog_name_galaxies= '{}{}'.format(FILE_NAME, type_galaxies)
        #terminal_galaxies= 'perl fiatfilter.pl -v "FWHM>{}*MAG_ISO+{} && FWHM>{} && CLASS_STAR<0.1 && FLUX_ISO<3000000" {}>{}'.format(m, n, FWHM_max_stars, catalog_name_good, catalog_name_galaxies)
        terminal_galaxies= 'perl fiatfilter.pl -v "FWHM>{}*MAG_ISO+{} && FWHM>{} && FLUX_ISO<3000000" {}>{}'.format(m, n, FWHM_max_stars, catalog_name_good, catalog_name_galaxies)
        subprocess.call(terminal_galaxies, shell=True)
        fcat_galaxies=np.genfromtxt(catalog_name_galaxies, names=names)
        #subprocess.call('./fiatreview {} {}'.format(fits, catalog_name_galaxies), shell=True)

        magnitude1='mag_iso'
        magnitude2='FWHM'
        plotter(fcat_good, magnitude1, magnitude2, 2, '$mag(iso)$', '$FWHM/pixels$')
        mag_th= np.linspace(1, 30, 1000)
        p = np.poly1d(fit)
        plt.plot(mag_th, p(mag_th), 'b-')
        plt.show()
        ellipticity(fcat_galaxies, 9)
        plt.show()

        #(9.2.): Checking GALAXIES CATALOG with Source Extractor Neural Network Output
        P.figure()
        P.hist(fcat_galaxies['class_star'], 50, normed=1, histtype='stepfilled')
        P.show(block=False)



        # (***) CHECKING FOR STARS // GALAXIES DIVISION

        weights_stars=np.ones_like(fcat_stars['class_star'])/len(fcat_stars['class_star'])
        weights_galaxies=np.ones_like(fcat_galaxies['class_star'])/len(fcat_galaxies['class_star'])
        weights_all = np.ones_like(fcat_good['class_star'])/len(fcat_good['class_star'])
        
        plt.figure()
        plt.hist(fcat_stars['class_star'], weights = weights_stars, bins= 5, histtype='stepfilled', label ='stars')
        plt.hist(fcat_galaxies['class_star'], weights = weights_galaxies, bins= 5, histtype='stepfilled', label ='galaxies')
        plt.legend(loc='upper right')
        plt.xlabel('$class_{star}$', labelpad=20, fontsize=20)
        plt.ylabel('$Frequency$', fontsize=20)
        plt.ylim(0,0.6)
        plt.show()
        plt.hist(fcat_good['class_star'], color= 'r', weights = weights_all, bins=50, histtype='stepfilled', label ='all')
        plt.legend(loc='upper right')
        plt.xlabel('$class_{star}$', labelpad=20, fontsize=20)
        plt.ylabel('$Frequency$', fontsize=20)
        plt.ylim(0,0.6)
        plt.show()

        plt.show()


        
        #(10): Calling Ellipto to recalculate shapes and ellipticities: ELLIPTO CATALOG
        print("")
        print("Now it is necessary to call ellipto in order to obtain in a proper way sizes and shapes both for galaxies and stars")
        print("")
        print("STARS")
        print("")
        catalog_name_ellipto_stars='{}{}'.format(FILE_NAME, type_ellipto_stars)
        ellipto_caller(catalog_name_stars, fits, catalog_name_ellipto_stars)
        print("GALAXIES")
        catalog_name_ellipto_galaxies='{}{}'.format(FILE_NAME, type_ellipto_galaxies)
        ellipto_caller(catalog_name_galaxies, fits, catalog_name_ellipto_galaxies)
        print("DONE")
        print("")
        
        #(11): Now we clasify the catalogs obtained with ellipto forcing fiat filter: SHAPES CATALOG
        print("Filtering good obtained celestial object from ellipto using fiatfilter...")
        print("")
        print("STARS")
        catalog_name_shapes_stars='{}{}'.format(FILE_NAME, type_shapes_stars)
        fiatfilter_errcode_stars='perl fiatfilter.pl -v "errcode<2" {}>{}'.format(catalog_name_ellipto_stars, catalog_name_shapes_stars)
        subprocess.call(fiatfilter_errcode_stars, shell=True)
        print("")
        print("GALAXIES")
        print("")
        catalog_name_shapes_galaxies='{}{}'.format(FILE_NAME, type_shapes_galaxies)
        fiatfilter_errcode_galaxies='perl fiatfilter.pl -v "errcode<2" {}>{}'.format(catalog_name_ellipto_galaxies, catalog_name_shapes_galaxies)
        subprocess.call(fiatfilter_errcode_galaxies, shell=True)
        print("DONE")
        print("")
        
        #(12): Recalculating ellipticities for stars
        print("I'm recalculating ellipticities of the new star set after being corrected by ellipto")
        names_ellipto = ["x", "y", "mag_iso", "median", "ixx", "iyy", "ixy", "a_input", "b_input", "theta", "ellipticity", "errcode", "sigsky", "size", "flux", "mean_rho_4th", "sigma_e", "wander"]
        fiat_shapes_stars= np.genfromtxt(catalog_name_shapes_stars, names=names_ellipto)
        ellipticity(fiat_shapes_stars, 15)
        plt.show()
        
        print "Show ellipticy as a function of x and y"
        plotter(fiat_shapes_stars, 'x', 'ellipticity', 2, '$x/pixels$', '$\epsilon$')
        plt.show()
        plotter(fiat_shapes_stars, 'y', 'ellipticity', 2, '$y/pixels$', '$\epsilon$')
        plt.show()

        fiat_shapes_galaxies= np.genfromtxt(catalog_name_shapes_galaxies, names=names_ellipto)
        ellipticity(fiat_shapes_galaxies, 15)
        plt.show(block=False)
        
        #(13): STARS--> you obtain two fitting both for ellip_1 and ellip_2
        print("")
        print("I'm performing a fitting of those ellipticities e_1 and e_2: both a simple 2D polynomial fitting and a 2D Legendre Polynomial fitting")
        print("")
        dlscombine_file_pol=''
        dlscombine_file_leg=''
        #Let's call the function fit_Polynomial from ellip_fitting3.py
        fitting_file_ellip_pol=ellip_fit.fit_Polynomial(FILE_NAME, fiat_shapes_stars)
        #Create file read by dlscombine
        dlscombine_file_pol=specfile(fits, fitting_file_ellip_pol, FILE_NAME)
        print("")
        #Let's call the function fit_Legendre from ellip_fitting3.py
        fitting_file_ellip_leg=ellip_fit.fit_Legendre(FILE_NAME, fiat_shapes_stars)
        #Create file read by dlscombine
        
        if filter=='r':
            dlscombine_file_leg=specfile_r(fits, fitting_file_ellip_leg, FILE_NAME)
        
        if filter=='z':
            dlscombine_file_leg=specfile_z(fits, fitting_file_ellip_leg, FILE_NAME)
        
        #(14): Let's call DLSCOMBINE to correct PSF anisotropies
        print("I'm correcting PSF anisotropies using dlscombine: BOTH FOR POL AND LEG FITTING")
        print("")
        fits_pol='{}_corrected_pol.fits'.format(FILE_NAME, FILE_NAME)
        dlscombine_call_pol='./dlscombine_pol {} {}'.format(dlscombine_file_pol, fits_pol)
        subprocess.call(dlscombine_call_pol, shell=True)
        fits_leg='{}_corrected_leg.fits'.format(FILE_NAME, FILE_NAME)
        dlscombine_call_leg='./dlscombine_leg {} {}'.format(dlscombine_file_leg, fits_leg)
        subprocess.call(dlscombine_call_leg, shell=True)
        
        #(15): Call again Source Extractor only for the Legendre Polynomial fitting
        print("I'm calling again SExtractor to obtain a new catalog from the corrected picture (only from the leg fitting)")
        print("")
        catalog_name_corrected=sex_caller_corrected(fits_leg, FILE_NAME)
    
        #(16): Transform .cat into .fcat (FIAT) for the corrected catalog
    
        catalog_name_fiat_corrected='{}_corrected.fcat'.format(FILE_NAME)
        transform_into_fiat_corrected='perl sex2fiat.pl {}>{}'.format(catalog_name_corrected, catalog_name_fiat_corrected)
        subprocess.call(transform_into_fiat_corrected, shell=True)
        print("")
        
        array_file_name.append(catalog_name_fiat_corrected)
        cont = cont + 1

    NAME_1= array_file_name[0]
    NAME_2= array_file_name[1]
    BEFORE_NAME_1 = NAME_1.find('.')
    FILE_NAME_1 = NAME_1[:BEFORE_NAME]
    BEFORE_NAME_2 = NAME_2.find('.')
    FILE_NAME_2 = NAME_2[:BEFORE_NAME]
    #CROSS-MATCHING
    catag_r = CatalogReader(array_file_name[0])
    catag_r.read()
    catag_z = CatalogReader(array_file_name[1])
    catag_z.read()
    crossmatching = CrossMatching(catag_r.fcat, catag_z.fcat)
    crossmatching.kdtree(n=1*1e-06)
    crossmatching.catalog_writter('2CM_{}'.format(FILE_NAME_1), compare = '1to2')
    print '\n'
    crossmatching.catalog_writter('2CM_{}'.format(FILE_NAME_2), compare = '2to1')
    
    FILE_NAME_FINAL = raw_input("Please, tell me the FINAL name: ")
    
    if crossmatching.cont1to2<crossmatching.cont2to1:
        catag_final_1 = CatalogReader('2CM_{}{}'.format(FILE_NAME_1, type_fcat))
        catag_final_1.read()
        catag_final_2 = CatalogReader('2CM_{}{}'.format(FILE_NAME_2, type_fcat))
        catag_final_2.read()
        crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
        crossmatching_final.kdtree(n=1*1e-06)
        crossmatching.catalog_writter('{}'.format(FILE_NAME_FINAL), compare = '1to2')

    if crossmatching.cont1to2>crossmatching.cont2to1:
        catag_final_1 = CatalogReader('2CM_{}{}'.format(FILE_NAME_1, type_fcat))
        catag_final_1.read()
        catag_final_2 = CatalogReader('2CM_{}{}'.format(FILE_NAME_2, type_fcat))
        catag_final_2.read()
        crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
        crossmatching_final.kdtree(n=1*1e-06)
        crossmatching.catalog_writter('{}'.format(FILE_NAME_FINAL), compare = '2to1')

    if crossmatching.cont1to2==crossmatching.cont2to1:
        catag_final_1 = CatalogReader('2CM_{}{}'.format(FILE_NAME_1, type_fcat))
        catag_final_1.read()
        catag_final_2 = CatalogReader('2CM_{}{}'.format(FILE_NAME_2, type_fcat))
        catag_final_2.read()
        crossmatching_final = CrossMatching(catag_final_1.fcat, catag_final_2.fcat)
        crossmatching_final.kdtree(n=1*1e-06)
        crossmatching.catalog_writter('{}'.format(FILE_NAME_FINAL), compare = '1to2')

    catalog_name_fiat_corrected_final = '{}{}'.format(FILE_NAME_FINAL, fcat)

    #(17): Transform again tshe corrected catalog into a GOOD catalog
    catalog_name_corrected_good= '{}{}'.format(FILE_NAME_FINAL, type_good)
    terminal_corrected_good= 'perl fiatfilter.pl "x>{} && x<{} && y>{} && y<{}" {}>{}'.format(xmin_good, xmax_good, ymin_good, ymax_good, catalog_name_fiat_corrected_final, catalog_name_corrected_good)
    subprocess.call(terminal_corrected_good, shell=True)

    FILE_NAME_CORRECTED='{}_corrected'.format(FILE_NAME_FINAL)

    #(18): STARS CATALOG again...
    print("Now we need to repeat the classification to obtain only galaxies and stars as we did before")
    print("")
    print("Let me show you again the FWHM vs MAG plot \n")
    print("")
    fcat_corrected=np.genfromtxt(catalog_name_corrected_good, names=names)
    plotter(fcat_corrected, 'mag_iso', 'FWHM', 3, '$mag(iso)$', '$FWHM$')
    plt.show(block=False)
    print("First stars...")
    print("")
    catalog_name_fiat_corrected_stars=''
    catalog_name_fiat_corrected_stars, FWHM_max_stars=stars_maker(catalog_name_corrected_good, FILE_NAME_CORRECTED)
    fcat_stars_corrected=np.genfromtxt(catalog_name_fiat_corrected_stars, names=names)
    ellipticity(fcat_stars_corrected, 20)
    plt.show(block=False)
    
    #(19): GALAXIES CATALOG again...
    print("")
    print("Second galaxies...")
    print("")
    catalog_name_fiat_corrected_galaxies=galaxies_maker(catalog_name_corrected_good, FILE_NAME_CORRECTED, FWHM_max_stars)
    fcat_galaxies_corrected=np.genfromtxt(catalog_name_fiat_corrected_galaxies, names=names)


    # (***) CHECKING FOR STARS // GALAXIES DIVISION

    weights_stars=np.ones_like(fcat_stars_corrected['class_star'])/len(fcat_stars_corrected['class_star'])
    weights_galaxies=np.ones_like(fcat_galaxies_corrected['class_star'])/len(fcat_galaxies_corrected['class_star'])
    weights_all = np.ones_like(fcat_corrected['class_star'])/len(fcat_corrected['class_star'])
        
    plt.figure()
    plt.hist(fcat_stars_corrected['class_star'], weights = weights_stars, bins= 10, histtype='stepfilled', label ='stars')
    plt.hist(fcat_galaxies_corrected['class_star'], weights = weights_galaxies, bins= 15, histtype='stepfilled', label ='galaxies')
    plt.legend(loc='upper right')
    plt.xlabel('$class_{star}$', labelpad=20, fontsize=20)
    plt.ylabel('$Frequency$', fontsize=20)
    plt.show()
    plt.hist(fcat_corrected['class_star'], color= 'r', weights = weights_all, bins=50, histtype='stepfilled', label ='all')
    plt.legend(loc='upper right')
    plt.xlabel('$class_{star}$', labelpad=20, fontsize=20)
    plt.ylabel('$Frequency$', fontsize=20)
    plt.show()


    #(20): ELLIPTO CATALOG and SHAPES CATALOG (only galaxies) again...
    catalog_name_ellipto_stars_corrected='{}{}'.format(FILE_NAME_CORRECTED, type_ellipto_stars)
    ellipto_caller(catalog_name_fiat_corrected_stars, fits, catalog_name_ellipto_stars_corrected)
    catalog_name_ellipto_galaxies_corrected='{}{}'.format(FILE_NAME_CORRECTED, type_ellipto_galaxies)
    ellipto_caller(catalog_name_fiat_corrected_galaxies, fits, catalog_name_ellipto_galaxies_corrected)
    catalog_name_shapes_galaxies_corrected='{}{}'.format(FILE_NAME_CORRECTED, type_shapes_galaxies)
    fiatfilter_errcode_galaxies_corrected='perl fiatfilter.pl -v "errcode<2" {}>{}'.format(catalog_name_ellipto_galaxies_corrected, catalog_name_shapes_galaxies_corrected)
    subprocess.call(fiatfilter_errcode_galaxies_corrected, shell=True)


    catalog_name_shapes_stars_corrected='{}{}'.format(FILE_NAME_CORRECTED, type_shapes_stars)
    fiatfilter_errcode_stars_corrected='perl fiatfilter.pl -v "errcode<2" {}>{}'.format(catalog_name_ellipto_stars_corrected, catalog_name_shapes_stars_corrected)
    subprocess.call(fiatfilter_errcode_stars_corrected, shell=True)


if __name__ == "__main__":
    main()
