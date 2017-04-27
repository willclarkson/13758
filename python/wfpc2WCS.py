#
# wcsWFPC2.py
#

# WIC 2017-04-27 - utilities to handle the world coordinate system
# info for an WFPC2 file downloaded from the STScI archive ("c0m.fits"
# format)

from astropy.io import fits
from astropy.wcs import WCS
import os, copy

import numpy as np

class WFPC2WCS(object):

    """Object holding the wfpc2 wcs and methods to transform positions.

    Each chip is given the WCS object using astropy methods."""

    def __init__(self, pathFits='', doTest=False, runOnInit=True, \
                 Verbose=False):

        self.pathFits = pathFits[:]

        # list of FITS headers
        self.LHdrs = []

        # list of WCS objects
        self.LWCS = []
        
        # control variables
        self.doTest = doTest
        self.Verbose = Verbose
        
        # Set the variables appropriate for testing if asked
        self.setupTest()

        if runOnInit:
            self.loadHeaders()
            self.wcsFromHeaders()
        
    def loadHeaders(self):

        """Imports the headers to the instance"""

        self.LHdrs= []

        if not os.access(self.pathFits, os.R_OK):
            if self.Verbose:
                print "wfp2WCS.importHeaders FATAL - cannot read path %s" \
                    % (self.pathFits)
            return

        LHDU = fits.open(self.pathFits)
        
        for iHdr in range(len(LHDU)):
            self.LHdrs.append(copy.deepcopy(LHDU[iHdr].header))

        LHDU.close()

    def wcsFromHeaders(self):

        """Having loaded the headers, convert them into WCS"""

        if len(self.LHdrs) < 1:
            return

        for iWCS in range(len(self.LHdrs)):
            try:
                thisWCS = WCS(self.LHdrs[iWCS])
            except:
                thisWCS = None

            # append this to the master WCS list
            self.LWCS.append(thisWCS)

    def doPix2World(self, xPix, yPix, iHDU=1):

        """Returns world coordinates for the chip of interest. We assume the
        pixel-coordinates have already been selected for chip number"""

        raRet = np.array([])
        deRet = np.copy(raRet)
        
        if iHDU > len(self.LWCS):
            if self.Verbose:
                print "WFPC2WCS.doPix2World WARN - iHDU > len(WCS): %i, %i" \
                    % (iHDU, len(self.LWCS))
            return raRet, deRet

        # otherwise, do the conversion
        raRet, deRet = self.LWCS[iHDU].all_pix2world(xPix, yPix, 0)

        return raRet, deRet

    def doPix2WorldMulti(self, xPix, yPix, extens, badVal=-99.9):

        """For multi-extension photometry. Given pixel positions and
        corresponding chip identifiers (counting from 1 ), this
        returns pixel positions translated to the sky for each chip.

        """

        raRet = np.zeros(np.size(xPix)) + badVal
        deRet = np.copy(raRet)

        # unique indices in the input
        lInds = np.unique(np.sort(extens))

        for iExt in lInds:
            bChip = extens == iExt
            if np.sum(bChip) < 1:
                continue

            # don't try to go beyond the length of the wcs list
            if iExt >= len(self.LWCS):
                continue

            # now compute
            thisWCS = self.LWCS[iExt]
            raRet[bChip], deRet[bChip] = thisWCS.all_pix2world(xPix[bChip], yPix[bChip], 0)

            if self.Verbose:
                print "Done iExt %i: %i objects" % (iExt, np.sum(bChip))

        return raRet, deRet
                
    def setupTest(self):

        """Set some defaults if testing"""

        if not self.doTest:
            return

        dirFits = '/home/wiclarks/Data/scratch/testBoresight/veronicaDOLPHOT/testGAIA'
        filFits = 'u49n1801r_c0m.fits'

        self.pathFits = '%s/%s' % (dirFits, filFits)

###

def TestLoad():

    """Test/debug routine"""

    WW = WFPC2WCS()
    WW.doTest = True
    WW.setupTest()

    WW.loadHeaders()

    WW.wcsFromHeaders()
    print len(WW.LWCS)
