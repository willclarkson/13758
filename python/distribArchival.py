#
# distribArchival.py
#

# WIC 2017-04-25 - distribute archival observations of a field in
# preparation for proper motion analysis

import os, sys, time, glob
from astropy.io import fits

def findFits(dirTop="", srchString='.fits'):

    if len(dirTop) < 1:
        dirTop = "/Volumes/Promise Pegasus/top/share/Data/HST/ngc6528_archival"


    # search for files
    lFits = glob.glob("%s/*%s*" % (dirTop, srchString) )

    return lFits

def getFileInfo(filePath=''):
    
    if len(filePath) < 1:
        return

    if not os.access(filePath, os.R_OK):
        return

    hdr = fits.getheader(filePath, 0)
    fSho = os.path.split(filePath)[-1]
    fStem = fSho.split('_')[0]
    
    try:
        propID = hdr['PROPOSID']
        instrume = hdr['INSTRUME']
        targname = hdr['TARGNAME'].replace('-','')
        exptime = str(hdr['EXPTIME']).replace('.','p')
    except:
        noID = True
        print fSho
        return

    print fSho, fStem, propID, instrume, targname, exptime
    
def go(srchStr='.fits'):

    """Wrapper - does the searching"""

    LFits = findFits(srchString=srchStr)

    print len(LFits)

    for iFil in range(3):
        getFileInfo(LFits[iFil])
