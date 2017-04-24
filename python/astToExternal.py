#
# astToExternal.py
#

# General-purpose module to match input photometry to astrometric
# catalog.

import os, time
from astropy.table import Table

import numpy as np

# scikits-learn
from sklearn.mixture import GMM

# for querying an external catalog
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

# for parsing box region for co-ordinate import
import astropy.units as u
import astropy.coordinates as coord

# for matching on the sphere
from scipy.spatial import cKDTree

from astropy.stats import sigma_clip

# for visualizing our dataset
import matplotlib.pylab as plt
plt.ion()

class ObsCat(object):

    """Object holding photometry catalog to match. Includes methods for
selection"""
    
    def __init__(self, \
                 pathPhot='ACS_WFPC2-2_F625W_chip1_PHOT_withFlags.fits', \
                 selectionLims={'goodPhot':[0,5],\
                                'Instrumental_VEGAMAG':[10., 23.]},\
                 Verbose=True, \
                 maxCatRows = 10000):

        # input file information
        self.pathPhot = pathPhot[:]

        # master-table of photometry
        self.tPhot = Table()
        self.bUse = np.array([])
        
        # column names for RA, DEC
        self.colRA = 'RA'
        self.colDE = 'DEC'

        # convenience variables
        self.ra = np.array([])
        self.de = np.array([])
        
        # center of the co-ordinate set
        self.cenRA = 0.
        self.cenDE = 0.

        # coordinate-aligned box widths
        self.boxWidth = 0.
        self.boxHeight = 0.

        self.boxUnits='deg'
        self.cooFrame='icrs'

        # list of external catalogs returned
        self.catalog='gaia'
        self.LExternalCats = []
        self.tAst = Table()

        # maximum number of rows to query
        self.maxCatRows = maxCatRows
        
        # column names in comparison catalog; convenience-views
        self.colRA_Ast = 'RAJ2000'
        self.colDE_Ast = 'DEJ2000'
        #self.colRA_Ast = 'RA_ICRS'
        #self.colDE_Ast = 'DE_ICRS'
        self.ra_Ast = np.array([])
        self.de_Ast = np.array([])
        
        # plottable bounding box
        self.plotBox = np.array([])
        
        # selection criteria
        self.selectionLims = selectionLims.copy()

        # Matching criteria - distance
        self.distMax3D = 0.1

        # For checking matches: mag columns
        self.colMag = 'Instrumental_VEGAMAG'
        self.colMag_Ast = 'f.mag'
        
        self.Verbose = Verbose
        
    def loadPhot(self):

        """Loads the photom file into memory."""

        if not os.access(self.pathPhot, os.R_OK):
            if self.Verbose:
                print "obsCat.loadPhot WARN - cannot load catalog %s" \
                    % (self.pathPhot)
            return

        # load the table, initialize the boolean
        self.tPhot = Table.read(self.pathPhot)
        self.bUse = np.repeat(True, len(self.tPhot))
            
    def selectByLimits(self):

        """Updates the selection boolean by the specified limits. Is
        conservative in the sense that if a keyword is missing from
        the table, a warning rather than exception is raised."""

        if np.sum(self.bUse) < 1:
            return

        for criterion in self.selectionLims.keys():
            limLo = self.selectionLims[criterion][0]
            limHi = self.selectionLims[criterion][1]

            if not criterion in self.tPhot.colnames:
                if self.Verbose:
                    print "obsCat.selectByLimits WARN - criterion %s not in the photometry table. Ignoring." % (criterion)
                continue

            self.bUse = (self.bUse) & \
                        (self.tPhot[criterion] >= limLo) & \
                        (self.tPhot[criterion] < limHi)

            self.setRADEC()
            
    def setRADEC(self):

        """Convenience method - sets the ra, dec"""

        if not self.colRA in self.tPhot.colnames:
            print "ObsCat.setRADEC WARN - RA column not present: %s" \
                % (self.colRA)
        else:
            self.ra = self.tPhot[self.colRA]

        if not self.colDE in self.tPhot.colnames:
            print "ObsCat.setRADEC WARN - DEC column not present: %s" \
                % (self.colDE)
        else:
            self.de = self.tPhot[self.colDE]

    def findAlignedBBox(self):

        """Determines the bounding rectangle aligned with the coordinates"""

        self.findFieldCenter()
        self.findBoundWidthHeight()
        self.repBBoxAsPoly()
        
    def findFieldCenter(self):

        """Finds the field center using the extrema of detected objects"""

        if np.size(self.ra) < 1 or np.size(self.de) < 1:
            if self.Verbose:
                print "ObsCat.findFieldCenter WARN - one of ra or de is missing."
            return

        self.cenRA = 0.5*(np.max(self.ra) + np.min(self.ra))
        self.cenDE = 0.5*(np.max(self.de) + np.min(self.de))
        
    def findBoundWidthHeight(self):

        """Determines the width and height of the bounding rectangle that is
        aligned with the co-ordinate frame (for external catalog
        queries)

        """

        self.boxWidth = np.max(self.ra) - np.min(self.ra)
        self.boxHeight = np.max(self.de) - np.min(self.de)

    def repBBoxAsPoly(self):

        """Convenience function - represents the bounding box as a set of
        coordinates for easy plotting with matplotlib.

        """
        
        vRA = self.cenRA + 0.5 * self.boxWidth \
              * np.array([-1., 1., 1., -1., -1.])
        vDE = self.cenDE + 0.5 * self.boxHeight \
              * np.array([-1., -1., 1., 1., -1.])

        self.plotBox = np.vstack(( vRA, vDE ))

    def queryExternal(self):

        """Runs the query on specified external catalog(s)."""

        cen = coord.SkyCoord(ra=self.cenRA, \
                             dec=self.cenDE, \
                             unit=(u.deg, u.deg), \
                             frame=self.cooFrame)

        wid = Angle(self.boxWidth, self.boxUnits)
        hgt = Angle(self.boxHeight, self.boxUnits)

        if self.Verbose:
            t0 = time.time()
            print "AstCat.doQuery INFO: Attempting query of %s..." \
                % (self.catalog)

        # Update the row limit for the catalog
        Vizier.ROW_LIMIT = int(self.maxCatRows)

        # WATCHOUT - the width and height seem to be flipped by
        # Vizier.query_region (using astroquery 0.3.4). The switching
        # of labels below is NOT a bug!!
        result = Vizier.query_region(cen, \
                                     width=hgt, \
                                     height=wid, \
                                     catalog=self.catalog)

        if self.Verbose:
            print "AstCat.doQuery INFO: found %i rows [0] in %.3f seconds" \
                % (len(result[0]), time.time() - t0)

        # add this to the list of catalogs
        self.LExternalCats = result
        self.tAst = self.LExternalCats[0]
        self.populateAstPosns()

    def populateAstPosns(self):

        """Populates the RA, DEC of the comparison astrometric catalog"""

        if not self.colRA_Ast in self.tAst.colnames:
            if self.Verbose:
                print "ObsCat.populateAstPosns WARN - RA not populated."
            return

        if not self.colDE_Ast in self.tAst.colnames:
            if self.Verbose:
                print "ObsCat.populateAstPosns WARN - DEC not populated."
            return

        self.ra_Ast = self.tAst[self.colRA_Ast]
        self.de_Ast = self.tAst[self.colDE_Ast]

    def equat2xyz(self, ra=np.array([]), de=np.array([]) ):

        """Utility - converts equatorial coordinates to xyz positions"""

        cosAlpha = np.cos(np.radians(ra))
        sinAlpha = np.sin(np.radians(ra))

        cosDelta = np.cos(np.radians(de))
        sinDelta = np.sin(np.radians(de))

        xyz = np.zeros( (np.size(cosAlpha), 3) )

        xyz[:,0] = cosDelta * cosAlpha
        xyz[:,1] = cosDelta * sinAlpha
        xyz[:,2] = sinDelta

        return xyz

    def matchOnSphere(self, nClip=5, minVarRatio=2.):

        """Matches the obs to the astrometric catalog"""

        xyzObs = self.equat2xyz(self.ra, self.de)
        xyzAst = self.equat2xyz(self.ra_Ast, self.de_Ast)

        bTry = np.copy(self.bUse)

        # Minimum matching distance
        distMax = np.copy(self.distMax3D)
        for iIter in range(nClip):
        
            kdt = cKDTree(xyzAst)
            dists, indices = kdt.query(xyzObs[bTry], \
                                       distance_upper_bound=distMax)

            dra = self.ra[bTry] - self.ra_Ast[indices]
            dde = self.de[bTry] - self.de_Ast[indices]
            dde *= np.cos(np.radians(self.de_Ast[indices]))

            # use mixture modeling on this space. Pre-convert degrees
            # to arcsec to get round precision limitations
            clf = GMM(2)
            delts = np.vstack(( np.asarray(dra)*3600., np.asarray(dde)*3600. ))
            delts = np.reshape(delts, (np.size(dra), 2) )

            clf.fit(delts)

            # pick the winner from the stddev
            devs = np.sqrt(np.sum(clf.covars_**2, axis=1))           

            if np.abs((devs[1] - devs[0])/np.min(devs)) < minVarRatio:
                print "Dropping out of iterations at %i - clumps similar" \
                    % (iIter)
                break
                
            labelWinner = np.argmin(devs)

            # predict membership
            labelPredicted = clf.predict(delts)
            gBroad = np.where(labelPredicted != labelWinner)[0]

            # update the boolean and try again
            gTry = np.where(bTry)[0]
            bTry[gBroad] = False

            continue

        # That clips it down... Now run again, this time on the
        # clipped informatoin (if over-clipped the population just
        # gets split in two...)

        kdt = cKDTree(xyzAst)
        dists, indices = kdt.query(xyzObs[bTry], \
                                   distance_upper_bound=distMax)

        dra = self.ra[bTry] - self.ra_Ast[indices]
        dde = self.de[bTry] - self.de_Ast[indices]
        dde *= np.cos(np.radians(self.de_Ast[indices]))

        delts = np.vstack(( np.asarray(dra)*3600., np.asarray(dde)*3600. ))
        delts = np.reshape(delts, (np.size(dra), 2) )
        labelPredicted = clf.predict(delts)
        
        # Find the median quantities, converted to WFC pixels
        dClip = sigma_clip(dists*206265., sigma=3., iters=5)
        print "Offsets, WFC pix: %.3f, %.3f" \
            % (np.median(dra[~dClip.mask])*3600. * 50., \
            np.median(dde[~dClip.mask])*3600. * 50. )

        figDelta = plt.figure(2)
        figDelta.clf()

        # it's useful to show the deltas
        ax2 = figDelta.add_subplot(222)
        blah2 = ax2.scatter(delts[:,0], delts[:,1], s=16, c=labelPredicted, \
                           edgecolor='none')
        figDelta.colorbar(blah2)
        ax2.set_title("Last iteration %i: Narrow=%i" % (iIter, labelWinner))
        
        ax3 = figDelta.add_subplot(223)
        #blah = ax3.scatter(delts[:,0], delts[:,1], s=16, c=labelPredicted, \
        #                   edgecolor='none')
        blah = ax3.scatter(dra*3600., dde*3600., s=16, \
                           edgecolor='None')
        ax3.set_xlabel('Offset, arcsec')
        ax3.set_ylabel('Offset, arcsec')
        ax3.set_title('Accepted points for offset')
        #figDelta.colorbar(blah)
        
        ax1 = figDelta.add_subplot(221, sharex=ax3, alpha=0.5)
        dum = ax1.hist(dra*3600., bins=50)

        ax4 = figDelta.add_subplot(224, sharey=ax3)
        dum = ax4.hist(dde*3600., bins=50, alpha=0.5, orientation='horizontal')

    def showSelectedPoints(self):

        """Debug routine - shows the points selected thus far"""

        if np.sum(self.bUse) < 1:
            return
        
        plt.figure(1)
        plt.clf()

        # convenience views
        ra = self.tPhot[self.colRA]
        de = self.tPhot[self.colDE]

        plt.scatter(ra[self.bUse], de[self.bUse])
        plt.xlabel(r"$\alpha$", fontsize=16.)
        plt.ylabel(r"$\delta$", fontsize=16.)

        # show comparison points
        if np.size(self.ra_Ast) > 0:
            plt.scatter(self.ra_Ast, self.de_Ast, c='r', marker='x')
        
        if np.size(self.plotBox) > 0:
            plt.plot(self.plotBox[0], self.plotBox[1], 'g-', lw=2)

def testRegion(magLo=14., magHi=18.):

    """Tests the region bounding the photometry catalog"""

    OC = ObsCat()
    OC.loadPhot()

    # adjust the selection limits if needed
    OC.selectionLims['Instrumental_VEGAMAG'] = [magLo, magHi]    
    OC.selectByLimits()
    
    OC.findAlignedBBox()
    OC.showSelectedPoints()

    OC.queryExternal()

    print OC.tAst.colnames
    
    #return
    
    print "Starting matching..."
    OC.matchOnSphere()
    
    OC.showSelectedPoints()

    
    return
    
    # Now build an astrometric catalog from the result and query it
    AC = AstCat(OC.cenRA, OC.cenDE, OC.boxWidth, OC.boxHeight)

    AC.doQuery()
    print AC.tCat.colnames
    
#    print np.shape(OC.bbox)
#    OC.findMBR()
#    print np.shape(OC.bbox)

#    print OC.bbox
    
#    OC.showSelectedPoints()
