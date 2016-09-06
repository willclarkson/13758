#
# dolXYtoRD.py
#

#
# WIC 2016-09-06
# 

# Purpose: translate pixel coordinates to world coordinates using FITS
# header

import os

from astropy.table import Table
from astropy.io import fits
from astropy import wcs

def xyToRd(filPhot='blah.xym', filWCS='', filOut='', \
               filInfo='blah.info', \
               colX='col3', colY='col4', \
               dirWCS='', Verbose=True, \
               Debug=True):

    """Converts pixels to ra, dec, using the World Coordinate System
    encoded in (usually fits) file filWCS. Returns the name of the output file if successful, a zero-length blank string otherwise.

    Inputs:

    filPhot -- Input photometry file with pixel coordinates
    
    filWCS -- Input file containing WCS header. 

    filOut -- Output file name to contain the X, Y, RA, DEC, Mag

    filInfo -- WCS file might not have a filename like filPhot. In
    this case, use the .info output from dolphot to find the name of
    the reference file.

    dirWCS -- Directory holding the WCS fits file.

    colX -- column in the photometry table holding the X pixel coords

    colY -- column in the photometry table holding the Y pixel coords

    Verbose -- Provide debug messages to the screen?

    Debug -- print info about the files without actually doing anything on the data

    --

    
    If filWCS has zero length, or is not found in the working
    directory, its name will be constructed assuming a convention."""

    ### NOTE - most of the lines in what follow are filename
    ### handling. The actual conversion doesn't start until the
    ### Table.read() syntax about halfway down.

    if not os.access(filPhot, os.R_OK):
        if Verbose:
            print "dolXYtoRD.xyToRd WARN - cannot read input path %s" \
                % (filPhot)
        return ''

    # Where are the WCS files to be found?
    if len(dirWCS) < 1:
        dirWCS = os.getcwd()
    else:
        if not os.access(dirWCS, os.R_OK):
            dirWCS = os.getcwd()
            if Verbose:
                print "dolXYtoRD.xyToRd WARN - WCS directory not readable: %s" \
                    % (dirWCS)
                print "dolXYtoRD.xyToRd WARN - Will search current directory for WCS"

    # Construct the WCS filename using our rules
    consWCS = getWcsNameFromPhot(filPhot)

    # If we've supplied an "info" file to get the reference image
    # name, use that instead of the photometry file name.
    if len(filInfo) > 0:
        if os.access(filInfo, os.R_OK):
            consWCS = getWcsNameFromInfo(filInfo)

    # If no WCS file was given, OR the given WCS file is not readable,
    # use the constructed WCS filename to search for the header.
    if len(filWCS) < 1:
        pathWCS = '%s/%s' % (dirWCS, consWCS)
    else:
        pathWCS = '%s/%s' % (dirWCS, filWCS)
        if not os.access(pathWCS, os.R_OK):
            pathWCS = '%s/%s' % (dirWCS, consWCS)

    # If we got here, then we have found both the photometry file and
    # the file holding the WCS information that we'll need. Construct
    # the output filename too.
    if len(filOut) < 4:  # (I like descriptive filenames)
        filOut = getOutNameFromPhot(filPhot)

    # some lines for debugging
    if Debug or Verbose:
        print "DBG: filPhot: ", filPhot, os.access(filPhot, os.R_OK)
        print "DBG: filInfo: ", filInfo, os.access(filInfo, os.R_OK)
        print "DBG: pathWCS: ", pathWCS, os.access(pathWCS, os.R_OK)
        print "DBG: filOut: ", filOut, os.access(filOut, os.R_OK)
        
        if Debug:
            return

    ## Having decided what to call everything, now look for them on
    ## disk

    # OK by this point we have our path for the WCS information. If
    # this isn't found on disk we should gracefully exit
    if not os.access(pathWCS, os.R_OK):
        if Verbose:
            print "dolXYtoRD.xyToRd FATAL - cannot read WCS path %s" % (pathWCS)
        return 


    # That's the hard part done. Now we do the easy part - read the
    # data, apply the transformation, and write the result.
    tPhot = Table.read(filPhot, format='ascii')
    
    hdulist = fits.open(pathWCS)

    #### NOTE we'll need to be a bit careful here... if using _flc
    #### files, then each chip will have a slightly different WCS
    #### header. So we'll have to find a way to account for the chips
    #### (should be easy). We need to know the convention!!
    w = wcs.WCS(hdulist[1].header, fobj=hdulist)  # need to get the extension right
    world = w.wcs_pix2world(tPhot[colX], tPhot[colY], 1)

    # close the file handle
    hdulist.close()

    tPhot['RA'] = world[0]
    tPhot['DEC'] = world[1]

    tPhot.write(filOut, overwrite=True)

    return filOut

### The following two methods construct the WCS filename and the
### output filename from the input filename. 

### You'll notice the two methods are basically identical, and indeed
### could have been written as a single method, just fed different
### arguments. However it's possible that the two cases might end up
### with different operations to construct the return values, so I
### decided to make them two separate routines for ease of reading and
### debugging!

def getWcsNameFromPhot(filPhot='blah.xym', sSrch='.xym', sRep='.fits'):

    """Turns the photometry name into the WCS header filename"""

    # NOTE - search and replace doesn't work if the output file has no
    # extension! Instead do this the long way...
    sStem = filPhot.split(sSrch)[0]
    outStr = '%s%s' % (sStem, sRep)

    
    return outStr

def getOutNameFromPhot(filPhot='blah.xym', sSrch='.xym', sRep='.ram'):

    """Turns the photometry name into the output filename"""

    # NOTE - search and replace doesn't work if the output file has no
    # extension! Instead do this the long way...
    sStem = filPhot.split(sSrch)[0]
    outStr = '%s%s' % (sStem, sRep)

    return outStr

def getWcsNameFromInfo(filInfo='blah.info', sSplit='.chip', sRep='_flc.fits'):

    """Gets the wcs file name from the .info file"""

    if not os.access(filInfo, os.R_OK):
        return 'NON'

    # Judging by the first example, the second line of the .info file
    # is the reference image.

    ### Change these rules to use a different strategy to find the
    ### reference image.
    with open(filInfo, 'r') as rObj:
        iCount=0
        for line in rObj:
            if iCount > 1:
                break
            sFil = line.strip()
            iCount = iCount + 1

    # Apply our rules to turn the first image name into the WCS image
    # name.
    fileStem = sFil.split(sSplit)[0]
    sWCS = '%s%s' % (fileStem, sRep)

    return sWCS
