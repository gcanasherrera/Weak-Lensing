# Weak-Lensing Package

Synopsis
====================
Weak Lensing is a repository that contains a set of programs plus the main script required to detect clusters of galaxies provided a .fits image of a sky region. The detection of these cluster of galaxies is performed thanks to the concept of Weak Lensing. Light coming from distant objects suffers a deviation due to the mass present in the universe.This deviation produces a shear in the size of distant objects which can be measured.

Documentation
====================
Weak-Lensing package is divided in different parts. The Main Script develops several steps required to measure the shear of ddistant objects. The procedure is based in 22 different steps, whose description is commenting in WL_Script.py. WL_Script has an auxiliary program that contains set of defined functions called by the Main Script (WL_Utils.py).Finally, the package contains three programs whose goals and relation with WL_Script.py are described below.

WL_Script
---------------------
### Dependences

### Description

#### Step 1: Define the ending of the input/output files

#### Step 2: Read the .fits image introduced in screen. Obtain the name of the .fits image, which is the base of the rest of the catalog and images' names.

#### Step 3: Call Source Extractor

#### Step 4: Transform the catalog provided by Source Extractor after analyzing the image (.cat) into a FIAT file (.fcat)

The FIAT File format 1.0 is defined [here](http://dls.physics.ucdavis.edu/fiat/fiat.html) by David Wittman. The transformation is performed thanks to a perl script converter called sex2fiat.pl, also included in the previous link.

#### Step 5: Read the new FIAT catalog
Read the FIAT files using the function genfromtxt in numpy library.

#### Step 6: Plot the size of the celestial objects as a function of their magnitude
Through this graph it is possible to stablish a classification of those celestial objects as galaxies or stars

#### Step 7: Obtain a catalog without blank spaces
When the x and y positions of the celestial objects extracted by Source Extractor are represented, it is possible to observe that there are some missing pixels in the image due to saturation of the camera. We avoid taking those bad pixels bonding the image. The result is a GOOD catalog with no blank spaces on it.

#### Step 8: 


> This is a blockquote.
> 
> This is the second paragraph in the blockquote.
>
> ## This is an H2 in a blockquote



