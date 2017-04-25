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

def getFileInfo(filePath='', doMove=False):
    
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
    except:
        noID = True
        print fSho
        return

    # if no exptime, must be calib
    try:
        exptime = str(hdr['EXPTIME']).replace('.','p')
        filtr = getFilter(hdr)
        isCalib = False
    except:
        isCalib = True

    if isCalib:
        filtr='CALIB'
        exptime = ''

    if len(str(propID)) < 1:
        propID = 'CAL'

    if len(targname) < 1:
        targname = 'NONE'
        
    #print fSho, fStem, propID, instrume, targname, exptime, filtr

    # use the information we have gathered to produce a destination
    # directory
    dirSub = '%s/%s/%s/%s' % (propID, targname, instrume, filtr)

    # drizzled files get their own directory - at least one
    # investigator here has combined multiple exptimes into the same
    # images.
    if fSho.find('drz.fits') < 0 and len(exptime) > 0:
        dirSub = '%s/%s' % (dirSub, exptime)
    else:
        dirSub = '%s/%s' % (dirSub, fSho.split('.fits')[0])

    print fSho, fStem, dirSub

    # ok now find all the fits files matching this particular stem
    dirIn = os.path.split(filePath)[0]
    dirOu = dirIn+'/'+dirSub

    LMove = glob.glob('%s/%s*fits*' % (dirIn, fStem))

    # now do the move
    for sFil in LMove:
        pathSrc = sFil[:]
        filDest = os.path.split(sFil)[-1]
        pathDest = '%s/%s' % (dirOu, filDest)
        # print pathSrc, pathDest

        if doMove:
            if not os.access(dirOu, os.R_OK):
                os.makedirs(dirOu)
            os.rename(pathSrc, pathDest)
        
def getFilter(hdr):

    """Returns the filter from the header"""

    filtr=''
    for filtKey in ['FILTER1', 'FILTNAM1', 'FILTER']:
        #if len(filtr) > 0:
        #    continue
        try:
            filtr=hdr[filtKey]
        except:
            notFilt = True

    if len(filtr) < 1:
        filtr='NONE'
            
    return filtr
            
def go(srchStr='.fits', doAll=True, doMove=False):

    """Wrapper - does the searching"""

    LFits = findFits(srchString=srchStr)

    print len(LFits)

    iMax = 3
    if doAll:
        iMax = len(LFits)
    
    for iFil in range(iMax):
        getFileInfo(LFits[iFil], doMove=doMove)

def goInSequence():

    """Go through the frames in sequence"""

    for sStr in ['drz.fits', 'flt.fits', 'c0m.fits']:
        go(sStr, True, True)
