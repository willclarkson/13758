## python ##

* dolXYtoRD.py -- Convert dolphot pixel coordinates to equatorial 
coordinates using astropy functionality.

* dolproc13758.py -- Converts dolphot text output into merged fits
  file for convenience, uses the WCS of a reference image to add RA,
  DEC columns.

Example callings (assuming the fits file below is the reference file
for the dolphot output given):

* cd /home/ferreir3/Desktop/PROJECT/NGC6553/Dolphot_Results/flc/F3336W

* import dolproc13758

* import astToExternal

* pathPhot = dolproc13758.loadAndProject('OUTPUT_WFC3_F336W_flc', 'icp003drq_flc.fits')

* astToExternal.testRegion(pathIn=pathPhot)
