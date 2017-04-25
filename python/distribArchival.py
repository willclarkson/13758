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

def go():

    """Wrapper - does the searching"""

    LFits = findFits()

    print len(LFits)
