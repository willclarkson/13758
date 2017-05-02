#
# dolproc.py
#

# WIC 2016-12-19

# convert raw DOLPHOT output into a form compatible with TOPCAT, etc.

from astropy.table import Table, Column
import os
import glob
import time

import copy

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

# convenience methods for handling multiple WCS for WFPC2
from wfpc2WCS import WFPC2WCS

class PhotPaths(object):

    """Utility object to find the DOLPHOT photometry output in a results
    directory. Also finds the parameters

    """
    
    def __init__(self, fieldSplit='09',\
                     cam='ACS',field='SWEEPS', filtr='F625W', \
                     Verbose=True, \
                     dirBase='/home/wiclarks/Data/HST/12020/howardba_results/results', setOnInit=True):

        self.dirBase=dirBase[:]

        self.cam = cam[:]
        self.field = field[:]
        self.filtr = filtr[:]
        self.fieldSplit=fieldSplit[:]
        self.dirPhotTop=''
        self.dirsPhot = []

        # initialise the top-level photometry directory and find the
        # photometry files
        if setOnInit:
            self.findPhotDirs()

        # file extensions we need in order to make progress
        self.extInfo = 'info'
        self.extCols = 'columns'

        # store the paths in a dictionary, with one keyword per chip
        self.pathsPhot={}
        
        self.Verbose=Verbose
        
    def setPhotTop(self):

        """Builds the photometry directory from the options chosen"""

        self.dirPhotTop='%s/%s/%s/%s' \
            % (self.dirBase, self.cam, self.field, self.filtr)

        if len(self.fieldSplit) > 1:
            self.dirPhotTop = '%s/%s' % (self.dirPhotTop, self.fieldSplit)

    def findPhotDirs(self):

        """Finds the photometry directories"""

        self.setPhotTop()
        self.dirsPhot = glob.glob('%s/chip*' % (self.dirPhotTop))

        #print self.dirPhotTop
        #print self.dirsPhot
        
    def populatePathsPhot(self):

        """Populates the paths for the photometry files"""

        for sChip in self.dirsPhot:
            filePhot, fileInfo, fileCols = self.getPhotFile(sChip)

            nFound = 0
            for sFil in [filePhot, fileInfo, fileCols]:
                if os.access(sFil, os.R_OK):
                    nFound = nFound + 1

                    
            # ignore this chip if <3 found
            if nFound < 3:
                continue

            # populate the paths dictionary. Use a shortened version
            # for the key, storing the fullpath inside it
            sKey = os.path.split(sChip)[-1]
            self.pathsPhot[sKey] = {'dirFull':sChip[:]}
            self.pathsPhot[sKey]['phot'] = os.path.split(filePhot)[-1]
            self.pathsPhot[sKey]['info'] = os.path.split(fileInfo)[-1]
            self.pathsPhot[sKey]['columns'] = os.path.split(fileCols)[-1]

            # add the path to the parameter file
            self.pathsPhot[sKey]['parsDir'] \
                = '%s/%s/%s/%s' \
                % (self.dirBase,self.cam, self.field, self.filtr)

            # predict the path to the drizzled file
            self.pathsPhot[sKey]['drzDir'] \
                = self.pathsPhot[sKey]['parsDir'][:].replace('/results', '/data')

    def getPhotFile(self, dirIn='', stem=''):

        """For a particular directory, finds the photometry file.

        User-input "stem" can be used to override the search result."""

        filePhot = ''
        fileInfo = ''
        fileCols = ''
        if not os.access(dirIn, os.R_OK):
            if self.Verbose:
                print "PhotPaths.getPhotFile WARN - input directory not readable: %s" % (dirIn)
            return filePhot, fileInfo, fileCols

        # find the stem if not already given
        if len(stem) < 1:
            LInfo = glob.glob('%s/*.%s' % (dirIn, self.extInfo))

            # return the initialized value if no info found
            if len(LInfo) < 1:
                return filePhot, fileInfo, fileCols

            # decide on the file stem to use. Slightly awkward syntax
            # here since we want to be able to specify which stem to
            # use as well.
            stem = os.path.splitext(os.path.split(LInfo[0])[-1])[0]

        # construct filenames from the stem (we'll search for them later)
        filePhot='%s/%s' % (dirIn, stem[:])
        fileInfo='%s.%s' % (filePhot, self.extInfo)
        fileCols='%s.%s' % (filePhot, self.extCols)

        return filePhot, fileInfo, fileCols

class Phot(object):

    """Class to read and translate DOLPHOT photometry output"""
    
    def __init__(self, dirPhot='', filPhot='', filInfo='', filCols='', \
                 dirRef='', filRef='', dirRefFits='', \
                 dirOut='', filOut='', \
                 filRefFits = 'NONE', \
                 Verbose=True):

        # location information for photometry
        self.dirPhot = dirPhot[:]
        self.filPhot = filPhot[:]
        self.filInfo = filInfo[:]
        self.filCols = filCols[:]

        # location information for the reference file
        self.dirRef = dirRef[:]
        self.filRef = filRef[:]
        self.dirRefFits = dirRefFits[:]

        # reference image containing the WCS of the reference image
        self.filRefFits = filRefFits[:]

        # output directory
        self.dirOut = dirOut[:]
        self.filOut = filOut[:]

        # flag - did we find the photometry file?
        self.foundPhot=False
        self.foundRef=False
        self.foundRefFits = False

        # is the reference fits image WFPC2?
        self.refIsWFPC2 = False
        
        # Column names for the photometry (since "col1, col2..."
        # aren't very readable)
        self.colNames = {}
        self.cols2Write = {}
        
        # table holding the photometry
        self.tPhot = Table()

        # WCS object
        self.wcs = None
        
        # control variable
        self.Verbose=Verbose
        
    def findPaths(self):

        """Checks that the needed photometry files are all present."""

        # carpentry... if the filPhot was given as a full path, not
        # just a directory...

        if self.filPhot.find('/') > -1:
            self.pathPhot = self.filPhot[:]
        else:
            self.pathPhot='%s/%s' % (self.dirPhot, self.filPhot)

        if self.filInfo.find('/') > -1:
            self.pathInfo = self.filInfo[:]
        else:
            self.pathInfo='%s/%s' % (self.dirPhot, self.filInfo)

        if self.filCols.find('/') > -1:
            self.pathCols = self.filCols[:]
        else:
            self.pathCols='%s/%s' % (self.dirPhot, self.filCols)

        nFound = 0
        #for sPath in [self.pathPhot, self.pathInfo, self.pathCols]:
        #    if os.access(sPath, os.R_OK):
        #        nFound = nFound + 1

        # allow updating by zipped file if present
        for sAttr in ['pathPhot', 'pathInfo', 'pathCols']:
            sPath = getattr(self, sAttr)
            if os.access(sPath, os.R_OK):
                nFound = nFound + 1
            else:
                sZip = '%s.gz' % (sPath)
                if os.access(sZip, os.R_OK):
                    setattr(self, sAttr, sZip)
                    nFound = nFound + 1
            
                
        if nFound > 2:
            self.foundPhot=True

        # now do the same with the reference file
        self.pathRef= '%s/%s' % (self.dirRef, self.filRef)
        if os.access(self.pathRef, os.R_OK):
            self.foundRef=True

        #print "INFO - pathRef: %s" % (self.pathRef)

    def parseParsFile(self, filExt='fits'):

        """Parses the parameter file to get the reference fits image used"""

        self.filRefFits = 'NONE'
        sSrch = 'img0_file'
        
        if not self.foundRef:
            return
        
        # now load the image
        with open(self.pathRef, 'r') as rObj:
            for sLine in rObj:
                if sLine.find(sSrch) < 0:
                    continue

                dum = sLine.strip().split(' ')[-1]
                self.filRefFits = '%s.%s' % (dum, filExt)

    def findRefImg(self):

        """Looks for the reference image """

        self.pathRefFits = False
        self.pathRefFits = '%s/%s' % (self.dirRefFits, self.filRefFits)

        if self.Verbose:
            print "dolproc.Phot.findRefImg INFO - looking for reffits: %s" \
                % (self.pathRefFits)

        # use glob to search for files with this or similar names
        if os.access(self.pathRefFits, os.R_OK):
            self.foundRefFits = True

    def loadPhot(self):

        """Loads the photometry file"""
                
        if not self.foundPhot:
            return

        self.setupPhotCols()
        
        if self.Verbose:
            print "dolproc.Phot.loadPhot INFO - loading dolphot file  %s ..." \
                % (self.pathPhot)
        t0 = time.time()
        self.tPhot = Table.read(self.pathPhot, format='ascii')
        if self.Verbose:
            print "dolproc.Phot.loadPhot INFO - ... done in %.1f seconds" \
                % (time.time() - t0)

        self.renamePhotCols()
        
    def setupPhotCols(self):

        """Helps set the column names for DOLPHOT output in a pithy but
human-readable way"""

        # I try to use the *.columns file as much as possible (in case
        # future versions reorder things for some reason).

        # Sometimes the second entry in the description is perfectly fine
        lJoin = [7, 8, 11, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        
        LCols = []
        LVals = []
        with open(self.pathCols, 'r') as rObj:
            for sRef in rObj:
                vThis = sRef.strip().split()

                # Sometimes the second descriptive word is the one we
                # want
                sCol = vThis[0].split('.')[0]
                iCol = int(sCol)

                # Which word do we lift?
                sVal = vThis[1][:]
                if iCol in lJoin:
                    sVal = '%s_%s' % (vThis[1], vThis[2])

                # remove any commas
                sVal = sVal.replace(',','')
                    
                LCols.append('col%i' % (iCol))
                LVals.append(sVal)

        # pass in the first 24 of these to the columns dictionary
        self.cols2Write = []
        for iCol in range(24):
            self.colNames[LCols[iCol]] = LVals[iCol]
            
        # now correct some of the less obvious ones by-hand:
        self.colNames['col3'] = 'X'
        self.colNames['col4'] = 'Y'
        self.colNames['col5'] = 'Chi_fit'
        self.colNames['col14'] = 'Normalized_rate'
        self.colNames['col15'] = 'Normalized_rate_uncertainty'

        # now set the columns to write
        self.cols2Write = []
        for iCol in range(1,25):
            sCol = 'col%i' % (iCol)
            self.cols2Write.append(self.colNames[sCol])
        
    def renamePhotCols(self):

        """Renames the photometry columns"""

        for sCol in self.colNames.keys():
            if not sCol in self.tPhot.colnames:
                continue

            self.tPhot[sCol].name = self.colNames[sCol]

    def setOutPath(self):

        """Sets the output path for fits"""

        # first off, don't let the output path be empty.
        if len(self.dirOut) < 3:
            self.dirOut = os.getcwd()

        # generate an output filename if none has been given
        if len(self.filOut) < 1:
            sSplit='/results/'
            sSplit = self.dirOut[:]+'/'
            sOut = self.dirPhot.split(sSplit)[-1].replace('/','_')
            self.filOut = '%s_PHOT.fits' % (sOut)

        # ensure the directory is meaningful...
        self.filOut = os.path.split(self.filOut)[-1]
            
        self.pathOut = '%s/%s' % (self.dirOut, self.filOut)

    def writePhot2Fits(self):

        """Writes the renamed photometry file to fits"""

        self.setOutPath()

        if self.Verbose:
            print "dolproc.Phot.writePhot2Fits INFO - writing to %s" \
                % (self.pathOut)
        
        # add the reference image as a keyword
        self.tPhot.meta['refPars'] = self.filRef[:]
        self.tPhot.meta['refImg'] = self.filRefFits[:]

        cols2Write = self.cols2Write[:]
        if len(cols2Write) < 1:
            cols2Write = self.tPhot.colnames
            
        self.tPhot[cols2Write].write(self.pathOut, format='fits', overwrite=True)

    def getPhotFromFits(self):

        """Used mainly in debugging the astrometry transformation. Reads the
photometry from fits file (since is about a factor 10 faster)

        """

        self.tPhot = Table.read('ACS_SWEEPS_F625W_09_chip1_PHOT.fits')        
        self.cols2Write = self.tPhot.colnames[0:24]
        

    def loadWCS(self):

        """Loads the WCS from the reference fits image"""

        if not self.foundRefFits:
            return

        self.checkIsWFPC2()

        if self.refIsWFPC2:
            if self.Verbose:
                print "dolproc.Phot.loadWCS INFO - assuming WFPC2 image."
            self.wcs = WFPC2WCS(self.pathRefFits, Verbose=False)
            return
        
        self.wcs = WCS(self.pathRefFits)

    def pix2Sky(self):

        """Creates standardized coords from the XRef and YRef values"""

        xPix = self.tPhot['X']
        yPix = self.tPhot['Y']

        if not self.refIsWFPC2:
            RA, DEC = self.wcs.all_pix2world(xPix, yPix, 0)
        else:
            # otherwise, loop through the chips
            xten = self.tPhot['Extension']
            RA, DEC = self.wcs.doPix2WorldMulti(xPix, yPix, xten)
            
            
        # just put in as a basic array for the moment... We'll worry
        # about units later if the calling program actually needs
        # them.
        self.tPhot['RA'] = RA
        self.tPhot['DEC'] = DEC

        # insert RA, DEC into the columns to write if not there already
        if not 'RA' in self.cols2Write:
            self.cols2Write.insert(4, 'RA')
        if not 'DEC' in self.cols2Write:
            self.cols2Write.insert(5, 'DEC')
            
    def checkIsWFPC2(self):

        """Checks if the reference image is WFPC2"""

        # we can update the method if we like here (e.g. check the
        # fits header for the INSTRUME keyword). For the moment, I
        # will just trust the STScI-set filename convention.

        self.refIsWFPC2 = False
        if os.path.split(self.pathRefFits)[-1].find('c0m') > -1:
            self.refIsWFPC2 = True 
        
def TestFindPhot(cam='ACS', field='SWEEPS', filtr='F625W', \
                     parsChip1='param1.pars', \
                 parsChip2='param09chip2.pars', \
                     fieldSplit='09', \
                     Debug=False, \
                     dirBase='/home/wiclarks/Data/HST/12020/howardba_results/results', \
                 fitsWithWCS=''):

    """Tests finding the photometry. Example call:

    

    If fitsWithWCS has nonzero length, doesn't use the paramfiles to
    search for the WCS, but instead takes in the fitsWithWCS
    file. WARNING - un-tested for flt.fits.

    """

    P = PhotPaths(cam=cam, field=field, filtr=filtr, fieldSplit=fieldSplit, \
                      dirBase=dirBase)
    P.populatePathsPhot()

    # look up the WCS fits file
    haveFitsWCS = False
    if len(fitsWithWCS) > 0:
        if os.access(fitsWithWCS, os.R_OK):
            haveFitsWCS = True
    
    # The parameter file names have probably been set by hand. Put
    # that into a dictionary we can conveniently use, later on.
    parFiles = {}
    parFiles['chip1'] =  parsChip1[:]
    parFiles['chip2'] =  parsChip2[:]
    
    for sChip in ['chip1', 'chip2']:

        paths = P.pathsPhot[sChip]
       
        PHOT = Phot(paths['dirFull'], paths['phot'], \
                        paths['info'], paths['columns'], \
                        paths['parsDir'])
        PHOT.filRef=parFiles[sChip]
        PHOT.dirRefFits = paths['drzDir'][:]
    
        # now see if we have everything
        PHOT.findPaths()  # find the photometry

        # only search for the reference fits file if asked. 
        if not haveFitsWCS:
            PHOT.parseParsFile()
            PHOT.findRefImg()
        else:
            PHOT.foundRefFits = True
            PHOT.pathRefFits = fitsWithWCS[:]

        # Generate a sensible output filename
        stemOut = '%s_%s_%s' % (P.cam, P.field, P.filtr)
        if len(P.fieldSplit) > 0:
            stemOut = '%s_%s' % (stemOut, P.fieldSplit)
        PHOT.filOut='%s_%s_PHOT.fits' % (stemOut, sChip)
        
        if Debug:
            PHOT.getPhotFromFits()
        else:
            PHOT.loadPhot()
            
        PHOT.loadWCS()
        PHOT.pix2Sky()
        PHOT.writePhot2Fits()

def loadAndProject(photStem='TEST', \
                       refFits='u49n1801r_c0m.fits', \
                       dirIn=os.getcwd(), \
                       dirRefFits=os.getcwd(), \
                       Verbose=True):

    """Loads dolphot photometry files beginning with "photStem,"
    merging into a .fits file for the photometry (currently in the
    frame of teh reference file only).

    Returns the output file written.

    If the reference fits image (or headerlet) with the world
    co-ordinate system is present (given by "refFits"), the world
    coordinates are added to the output fits file. Arguments:

    photStem -- filename stem for the DOLPHOT output

    refFits -- fits file with WCS info

    dirIn -- directory holding the photometry

    dirRefFits -- directory holding the reference fits file

    """

    # find the files
    PP = PhotPaths(setOnInit=False)
    fPhot, fInfo, fCols = PP.getPhotFile(dirIn, photStem)

    PHO = Phot(dirIn, fPhot, fInfo, fCols, \
               dirRefFits=dirIn, \
               filRefFits=refFits[:], \
               dirRef=dirRefFits, \
               filOut='%s_PHOT.fits' % (fPhot), \
               Verbose=Verbose)

    # Now try loading
    PHO.findPaths()
    PHO.loadPhot()

    # try loading the WCS and projecting coords onto the sky
    PHO.findRefImg()
    PHO.loadWCS()
    try:
        PHO.pix2Sky()
    except:
        print "LoadPhotom WARN - world co-ordinate transformation failed."

    # now write the output to fits
    PHO.writePhot2Fits()
    
    # try writing to text for easy RA, DEC plotting with Starlink-GAIA
    # PHO.tPhot[PHO.cols2Write].write('%s_PHOT.txt' % (photStem), format='ascii')

    return PHO.pathOut
