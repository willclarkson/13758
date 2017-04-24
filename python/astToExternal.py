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

# fits imput
from astropy.io import fits

# for visualizing our dataset
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
plt.style.use('ggplot')
plt.ion()

class ObsCat(object):

    """Object holding photometry catalog to match. Includes methods for
selection"""
    
    def __init__(self, \
                 pathPhot='ACS_WFPC2-2_F625W_chip2_PHOT_withFlags.fits', \
                 selectionLims={'goodPhot':[0,5],\
                                'Instrumental_VEGAMAG':[10., 23.]},\
                 Verbose=True, \
                 maxCatRows = 50000, \
                 widthMax=0.1, heightMax=0.1, \
                 ForceQuery=False):

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

        # min-max box widths
        self.widthMax = widthMax
        self.heightMax = heightMax

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

        # File to store external catalog
        self.dirAst='./tmpAST'
        self.filAst = 'TEST.fits'
        self.setPathAst()

        # Force re-query even if the file is already present?
        self.ForceQuery = ForceQuery
        
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
        self.distMax3D = 0.01

        # hard limits on deltas
        self.maxDeltaArcsec = 999.

        # Min needed to try matching
        self.minForMatch = 25

        # Number of passes at trimming when matching
        self.matchCleanIters = 5

        # For checking matches: mag columns
        self.colMag = 'Instrumental_VEGAMAG'
        self.colMag_Ast = 'magG'

        # some program flow 
        self.trimAst = False
        self.Verbose = Verbose
       
        # test nudge for debug
        self.nudgeRAtest = 0.5 / 3600.
        self.nudgeDEtest = -1.5 / 3600.

        # Need to do the offsets here...
        self.deltaRA = 0.
        self.deltaDE = 0.

        # Indices for matched objects in the photom and astrom
        # catalogs
        self.rowsPhot = np.array([])
        self.rowsAst = np.array([])
        
        # crval for deprojection
        self.crval = np.zeros(2)

        # columns for the focal plane coordinates
        self.colEpsilon = 'EPSILON'
        self.colEta = 'ETA'

    def loadPhot(self):

        """Loads the photom file into memory."""

        if not os.access(self.pathPhot, os.R_OK):
            if self.Verbose:
                print "obsCat.loadPhot WARN - cannot load catalog %s" \
                    % (self.pathPhot)
            return

        # load the table, initialize the boolean
        try:
            self.tPhot = Table.read(self.pathPhot)
        except:
            self.tPhot = Table.read(self.pathPhot, format='ascii')
        self.bUse = np.repeat(True, len(self.tPhot))

        # reference RA, DEC for tangent plane projection. We'll work
        # out how to populate this later.
        self.crvalInit = np.array([0., 0.])
        self.crvalCorr = np.array([0., 0.])

        # tangent-plane coordinates
        self.tpPhot = Table()
        self.tpAst = Table()
        
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

            if self.Verbose:
                print "obsCat.selectByLimits INFO - using criterion %.2f <= %s < %.2f" % (limLo, criterion, limHi)
            
            self.bUse = (self.bUse) & \
                        (self.tPhot[criterion] >= limLo) & \
                        (self.tPhot[criterion] < limHi)

        # whoops - bad indent
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

    def nudgeObsPositions(self):

        """Useful debug check: nudge the positions by some delta and re-fit"""
        
        self.ra += np.repeat(self.nudgeRAtest, np.size(self.ra))
        self.de += np.repeat(self.nudgeDEtest, np.size(self.ra))

        self.tPhot[self.colRA] += self.nudgeRAtest
        self.tPhot[self.colDE] += self.nudgeDEtest

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

        print "DBG:", np.shape(self.ra)
        
        self.boxWidth = np.max(self.ra) - np.min(self.ra)
        self.boxHeight = np.max(self.de) - np.min(self.de)

        # set width limits if we have them
        self.boxWidth = np.min([self.boxWidth, self.widthMax])
        self.boxHeight = np.min([self.boxHeight, self.heightMax])

    def repBBoxAsPoly(self):

        """Convenience function - represents the bounding box as a set of
        coordinates for easy plotting with matplotlib.

        """
        
        vRA = self.cenRA + 0.5 * self.boxWidth \
              * np.array([-1., 1., 1., -1., -1.])
        vDE = self.cenDE + 0.5 * self.boxHeight \
              * np.array([-1., -1., 1., 1., -1.])

        self.plotBox = np.vstack(( vRA, vDE ))

    def getAstCat(self, catName=''):

        """Obtains the astrometric catalog for the photfile of
interest. Optional input argument is the path to the astFile (if we
already have a file in mind, say).

        """

        # set the astrometric path to defaults
        if len(catName) < 1:
            self.setPathAst()
        else:
            self.pathAst = catName[:]
        
        if os.access(self.pathAst, os.R_OK) and not self.ForceQuery:

            if self.Verbose:
                print "getAstCat INFO - found previous ast cat %s. Loading." \
                    % (self.pathAst) 

            try:
                self.tAst = Table.read(self.pathAst)
            except:
                self.tAst = Table.read(self.pathAst, format='ascii.csv')

            
        if len(self.tAst.colnames) < 1:
            self.queryExternal()

        if self.trimAst:
            self.trimExternalCat(0.1)

        self.populateAstPosns()
            
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
            print "AstCat.doQuery INFO: Width %.2f', height %.2f'" \
                % (self.boxWidth*60., self.boxHeight*60.)
            
            print "AstCat.doQuery INFO: querying '%s' catalog..." \
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
            print "AstCat.doQuery INFO: found %i rows [0] in %.2f seconds" \
                % (len(result[0]), time.time() - t0)

        # add this to the list of catalogs
        self.LExternalCats = result
        self.tAst = self.LExternalCats[0]

        # these have been moved up to getAstCat
        #if self.trimAst:
        #    self.trimExternalCat(0.1)

        #self.populateAstPosns()

    def trimExternalCat(self, fBri=0.1):

        """Trim external catalog by some criterion"""

        if fBri >= 1.0:
            return

        if not self.colMag_Ast in self.tAst.colnames:
            if self.Verbose:
                print "ObsCat.trimExternalCat WARN - brightness column %s not found in astrometric catalog."
            return
            
        nOrig = len(self.tAst)

        vMag = self.tAst[self.colMag_Ast]
        lTrim = np.argsort(vMag)[0:np.int(nOrig*fBri)]
        self.tAst = self.tAst[lTrim]

        if self.Verbose:
            print "ObsCat.trimExternalCat INFO - trimmed Astr from %i to %i" \
                % (nOrig, len(self.tAst))

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

    def matchOnSphere(self, nClip=5, minVarRatio=2., showGroups=False):

        """Matches the obs to the astrometric catalog"""

        xyzObs = self.equat2xyz(self.ra, self.de)
        xyzAst = self.equat2xyz(self.ra_Ast, self.de_Ast)

        bTry = np.copy(self.bUse)

        # Minimum matching distance
        distMax = np.copy(self.distMax3D)
        for iIter in range(self.matchCleanIters):
        
            kdt = cKDTree(xyzAst)
            dists, indices = kdt.query(xyzObs[bTry], \
                                       distance_upper_bound=distMax)

            # print np.max(dists)*206265.

            # WARN - somehow, for Chandra, max(indices) =
            # np.shape(self.ra_Ast), which should be illegal. For the
            # moment, trim straight away. This might be a bug in kdtree?
            gBadIndices = np.where(indices >= np.size(self.ra_Ast))[0]
            indices[gBadIndices] = 0

            dra = self.ra[bTry] - self.ra_Ast[indices]
            dde = self.de[bTry] - self.de_Ast[indices]
            dde *= np.cos(np.radians(self.de_Ast[indices]))

            if np.abs(self.maxDeltaArcsec) < 100.:

                #print "DBG: ", iIter, np.sum(bTry), np.shape(dra), np.shape(bTry)
                gBad = np.where((np.abs(dra)*3600. > self.maxDeltaArcsec) | (np.abs(dde) * 3600. > self.maxDeltaArcsec))[0]

                bTry[gBad] = False

                # print "DBG:", iIter, np.sum(bTry)

                dists, indices = kdt.query(xyzObs[bTry], \
                                           distance_upper_bound=distMax)

                dra = self.ra[bTry] - self.ra_Ast[indices]
                dde = self.de[bTry] - self.de_Ast[indices]
                dde *= np.cos(np.radians(self.de_Ast[indices]))

            # break out if we've chopped down too many objects
            if np.sum(bTry) < self.minForMatch:
                break

            # use mixture modeling on this space. Pre-convert degrees
            # to arcsec to get round precision limitations
            clf = GMM(2)
            delts = np.vstack(( np.asarray(dra)*3600., np.asarray(dde)*3600. ))
            delts = np.reshape(delts, (np.size(dra), 2) )

            clf.fit(delts)

            # pick the winner from the stddev
            devs = np.sqrt(np.sum(clf.covars_**2, axis=1))           
            labelWinner = np.argmin(devs)
            
            if np.abs((devs[1] - devs[0])/np.min(devs)) < minVarRatio:
                print "Dropping out of iterations at %i - clumps similar" \
                    % (iIter)
                print clf.covars_
                print clf.means_
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

        # Since clipping is done, we can now pass these up to the
        # instance.
        self.rowsPhot = np.where(bTry)[0]
        self.rowsAst = np.copy(indices)
        
        # SEPARATE tangent-plane variables to be used when fitting
        # various things out...
        dde_tan = dde * np.cos(np.radians(self.de_Ast[indices]))
        
        delts = np.vstack(( np.asarray(dra)*3600., np.asarray(dde_tan)*3600. ))
        delts = np.reshape(delts, (np.size(dra), 2) )
        labelPredicted = clf.predict(delts)

        # Find the marginal values
        meanRA, meanDE = self.fitClumpFromMarginals(dra*3600., dde*3600., 3)

        print "Offsets, arcsec: %.3f, %.3f" % (meanRA, meanDE)

        # Find the median quantities, converted to WFC pixels
        dClip = sigma_clip(dists*206265., sigma=3., iters=5)
        print "Clipped medians, arcsec: %.3f, %.3f" \
            % (np.median(dra[~dClip.mask])*3600., \
            np.median(dde[~dClip.mask])*3600. )

        #self.deltaRA = np.median(dra[~dClip.mask])
        #self.deltaDE = np.median(dde[~dClip.mask])

        self.deltaRA = meanRA / 3600.
        self.deltaDE = meanDE / 3600.

        print "=== INFO === sample for fitting:", \
            np.shape(dists), np.sum(bTry), np.shape(self.rowsPhot), np.shape(self.rowsAst)

        ############ Deubg plot figure comes below. This is all very
        ############ messy still!!

        coloCrosshair='0.6'
        coloFound='m'
        lwFound=2
        nBins=150
        
        figDelta = plt.figure(2)
        figDelta.clf()

        figDelta.subplots_adjust(hspace=0., wspace=0.)

        # it's useful to show the deltas
        if showGroups:
            ax2 = figDelta.add_subplot(222)
            blah2 = ax2.scatter(delts[:,0], delts[:,1], s=9, c=labelPredicted, \
                                    edgecolor='none')
        #figDelta.colorbar(blah2)
        #ax2.set_title("Last iteration %i: Narrow=%i" % (iIter, labelWinner))
        
        ax3 = figDelta.add_subplot(223)
        #blah = ax3.scatter(delts[:,0], delts[:,1], s=16, c=labelPredicted, \
        #                   edgecolor='none')
        blah = ax3.scatter(dra*3600., dde*3600., s=3, \
                           edgecolor='None', alpha=1., color='b', zorder=10, \
                               norm=LogNorm())
        ax3.set_xlabel(r"$\Delta \alpha,  arcsec$")
        ax3.set_ylabel(r"$\Delta \delta,   arcsec$")
        # ax3.set_title('Accepted points for offset')

        #figDelta.colorbar(blah)
        ax3.set_xlim(-2., 2.)
        ax3.set_ylim(-2., 2.)
        ax3.grid(which='both', color='0.75')

        # plot guidelines for the axis
        ax3.plot([0.,0.], [-2., 2.], color=coloCrosshair, \
                 zorder=1)
        ax3.plot([-2., 2.], [0., 0.], color=coloCrosshair, \
                 zorder=1)

        # also plot the deltas found.
        ax3.plot([meanRA, meanRA], [-2., 2.], zorder=15, lw=lwFound, \
                 color=coloFound)
        ax3.plot([-2., 2.], [meanDE, meanDE], zorder=15, lw=lwFound, \
                 color=coloFound)

        #hist2 = ax3.hist2d(dra*3600., dde*3600., bins=(20,20), zorder=1, \
#range=[[-2., 2.],[-2., 2.]], norm=LogNorm()re)

        # Find the marginal values
        meanRA, meanDE = self.fitClumpFromMarginals(dra*3600., dde*3600., 2)

        print meanRA, meanDE

        ax1 = figDelta.add_subplot(221, sharex=ax3, alpha=0.5)
        dum = ax1.hist(dra*3600., bins=nBins, alpha=0.5, log=True, \
                           edgecolor='None', color='b')
        ax1.grid(which='both', color='0.75')
        ax1.xaxis.tick_top()

        # get the current limit
        ymax1 = ax1.get_ylim()[-1]
        
        ax1.plot([0.0,0.0], [1, ymax1]  , color=coloCrosshair)
        ax1.plot([meanRA, meanRA], [1, ymax1], lw=lwFound, \
                 color=coloFound)
        
        ax1.set_ylabel(r"$N(\Delta \alpha)$")


        ax4 = figDelta.add_subplot(224, sharey=ax3)
        dum = ax4.hist(dde*3600., bins=nBins, alpha=0.5, \
                           orientation='horizontal', log=True, \
                           edgecolor='None', color='b')
        ax4.grid(which='both', color='0.75')
        ax4.yaxis.tick_right()

        ax4.set_xlabel(r"$N(\Delta \delta)$")

        # get the current limit
        xmax4 = ax4.get_xlim()[-1]
        
        ax4.plot([1, xmax4], [0., 0.], color=coloCrosshair)
        ax4.plot([1, xmax4], [meanDE, meanDE]  , lw=lwFound, \
                 color=coloFound)


        ax3.set_xlim(-2., 2.)
        ax3.set_ylim(-2., 2.)

        figDelta.savefig('ast2ext_deltas.png')

    def fitClumpFromMarginals(self, \
                              dra=np.array([]), \
                              dde=np.array([]), \
                              nComps=3):

        """Fit the marginal distributions and find the clump where they
        cross

        """

        # Write out the pieces... we can listify later
        deltRA = np.reshape(dra, (np.size(dra), 1))
        deltDE = np.reshape(dde, (np.size(dde), 1))

        clfRA = GMM(nComps)
        clfDE = GMM(nComps)
        
        clfRA.fit(deltRA)
        clfDE.fit(deltDE)

        # Translate into 1D for convenience
        covRA = clfRA.covars_[:,0]
        covDE = clfDE.covars_[:,0]
        meanRA = clfRA.means_[:,0]
        meanDE = clfDE.means_[:,0]

        # Sort by increasing variances, identify the minimum-variance
        # component
        iRA = np.argsort(covRA)[0]
        iDE = np.argsort(covDE)[0]

        return meanRA[iRA], meanDE[iDE]
    
    def showDistributions(self, doLog=False):

        """Shows the distributions of all and matched objects"""

        if len(self.rowsPhot) < 1:
            return

        plt.figure(3, figsize=(8,8))
        plt.clf()
        
        print self.tPhot.colnames
        print len(self.rowsPhot), len(self.tPhot)

        mags = self.tPhot[self.colMag]
        bGood = mags < 35.

        rang = [np.min(mags[bGood]), np.max(mags[bGood])]

        dum1 = plt.hist(mags, bins=100, range=rang, color='0.8', \
                            log=doLog, edgecolor='0.7')
        dum1b = plt.hist(mags, bins=100, range=rang, color='k', \
                            log=doLog, fill=False, \
                             histtype='step')

        dum2 = plt.hist(mags[self.rowsPhot], bins=100, range=rang, \
                            log=doLog, color='b', edgecolor='0.7')

        plt.grid(which='both', color='0.75')

        plt.xlabel(self.colMag)
        plt.ylabel('N(%s)' % (self.colMag))

        # save to figure
        fignam='ast2ext_cdf_lin.png'
        if doLog:
            fignam = 'ast2ext_cdf_log.png'

        plt.savefig(fignam)

    def ensureCRVALset(self, Clobber=False, tol=1e-8):

        """Ensures we have something meaningful for the CRVAL."""

        # self.cenRA, self.cenDE
        if np.sum(np.abs(self.crval)) < tol and not Clobber:
            return

        if np.sum(np.abs([self.cenRA, self.cenDE])) < tol:
            self.findFieldCenter()

        self.crval[0] = np.copy(self.cenRA)
        self.crval[1] = np.copy(self.cenDE)

    def celestialToTP(self, ra, dec, ra0, de0):

        """Projects celestial coordinates onto the tangent plane"""

        self.ensureCRVALset()

        # right ascension only appears as deltas
        dra_rad = np.radians(ra - ra0)
        cosDe0 = np.cos(np.radians(de0))
        sinDe0 = np.sin(np.radians(de0))

        sinDe = np.sin(np.radians(dec))
        cosDe = np.cos(np.radians(dec))

        # denominator...
        denom = sinDe0 * sinDe + cosDe0*cosDe * np.cos(dra_rad)

        epsilon = (cosDe * np.sin(dra_rad)) / denom
        eta = (cosDe0*sinDe - sinDe0*cosDe*np.cos(dra_rad) ) / denom

        return np.degrees(epsilon), np.degrees(eta)

    def tableToTP(self, phot=True):

        """Translates a table to tangent plane coords."""
        
        colRA = self.colRA[:]
        colDE = self.colDE[:]
        if phot:
            tData = self.tPhot
        else:
            tData = self.tAst
            colRA = self.colRA_Ast[:]
            colDE = self.colDE_Ast[:]

        if not colRA in tData.colnames:
            return

        ra = tData[colRA]
        de = tData[colDE]
        epsilon, eta = self.celestialToTP(ra, de, self.crval[0], self.crval[1])
        
        tData[self.colEpsilon] = np.copy(epsilon)
        tData[self.colEta] = np.copy(eta)

    def deprojectBothToTangentPlane(self):

        """Convenience routine - deprojects both tables to the tangent
        plane"""

        # Should come up with a more pithy name for this method

        # ensure CRVAL is set
        self.ensureCRVALset()

        self.tableToTP(phot=True)
        self.tableToTP(phot=False)

    def showSelectedPoints(self):

        """Debug routine - shows the points selected thus far"""

        if np.sum(self.bUse) < 1:
            return
        
        plt.figure(1)
        plt.clf()

        # convenience views
        ra = self.tPhot[self.colRA]
        de = self.tPhot[self.colDE]

        plt.scatter(ra[self.bUse], de[self.bUse], s=4)
        plt.xlabel(r"$\alpha$", fontsize=16.)
        plt.ylabel(r"$\delta$", fontsize=16.)

        # show comparison points
        if np.size(self.ra_Ast) > 0:
            plt.scatter(self.ra_Ast, self.de_Ast, c='r', marker='x', s=4)
        
        if np.size(self.plotBox) > 0:
            plt.plot(self.plotBox[0], self.plotBox[1], 'g-', lw=2)

    def setPathAst(self):

        """Sets the path for the astrometric catalog"""

        if len(self.dirAst) < 1:
            self.dirAst = os.getcwd()

        filPhot = os.path.split(self.pathPhot)[-1]
        stemPhot = os.path.splitext(filPhot)[0]

        # add the selection criteria to the filename
        try:
            keys = self.selectionLims.keys()
            for sKey in keys:
                dThis = self.selectionLims[sKey]
                sSel = '%s_%.1f-%.1f' % (sKey, dThis[0], dThis[1])

                stemPhot = '%s_%s_%s' % (stemPhot, self.catalog, sSel)
        except:
            dumdum = 1
                
        self.pathAst = '%s/%s_AST.fits' % (self.dirAst, stemPhot)
            
    def writeAstCat(self, Clobber=False):

        """Writes astrom catalog to disk"""

        if len(self.tAst) < 1:
            return

        if os.access(self.pathAst, os.R_OK) and not Clobber:
            return
        
        # ensure the output directory exists. Inherit from the path
        # rather than the separate directory variable
        dirAst = os.path.split(self.pathAst)[0]
        if not os.access(dirAst, os.R_OK) and len(dirAst) > 1:
            os.makedirs(dirAst)

        if self.pathAst.find('csv') > -1:
            print "Writing astrom catalog..."
            
            self.tAst['RAJ2000', 'DEJ2000', 'magRF'].write(self.pathAst, format='ascii.csv')
            return
            
        self.tAst.write(self.pathAst, overwrite=True)
        
    def updateImageRefCoo(self, pathImg=''):

        """Utility routine - updates the header of a given image following the
        offsets found thus far.

        """

        if len(pathImg) < 1:
            return

        if not os.access(pathImg, os.R_OK):
            return

        # construct the output image
        pathOut = '%s_tweaked.fits' % (os.path.splitext(pathImg)[0])

        hdulist = fits.open(pathImg)
        
        hdulist[0].header['CRVAL1'] -= self.deltaRA
        hdulist[0].header['CRVAL2'] -= self.deltaDE

        hdulist.writeto(pathOut, clobber=True)
        hdulist.close()
        
        
def testRegion(magLo=14., magHi=18., pathIn='', tryChandra=False, \
               debugNudge=False, iters=5, colRA='RA', colDE='DEC', \
               colMag='Instrumental_VEGAMAG', \
               img2Shift='', showLog=False):

    """Tests the region bounding the photometry catalog.

    Example call for checking against sextractor magnitudes:


    astToExternal.testRegion(pathIn='f625w_detection.FIT', colRA='ALPHA_J2000', colDE='DELTA_J2000', colMag='MAG_ISO', magLo=-11., magHi=-8., iters=7, img2Shift='f625w_drz_sci_0p7.fits')
    """

    # note (to turn into proper documentation later) - to run this on
    # F390M testing on desktop:
    #
    # cd /home/wiclarks/Data/scratch/testBoresight/veronicaDOLPHOT

    # astToExternal.testRegion(pathIn='TEST_bothChips.fits', magHi=23, magLo=18, showLog=True)

    OC = ObsCat(pathPhot=pathIn[:])
    OC.colRA = colRA[:]
    OC.colDE = colDE[:]
    
    OC.matchCleanIters = iters

    OC.nudgeRAtest = 0.
    OC.nudgeDEtest = 0.
    
    if tryChandra:
        OC.pathPhot = 'psf_95p_p1.5keV.ph'
        OC.colMag = 'SRCH'
        OC.widthMax = 0.08
        OC.heightMax = 0.08
        magLo = -999.
        magHi =  999.
        OC.maxCatRows = 50000
        OC.trimAst = False
        OC.distMax3D = 0.1 #250. / 206265.
        OC.maxDeltaArcsec = 10.
        OC.matchCleanIters = 5
    OC.loadPhot()
    
    if debugNudge:
        OC.nudgeObsPositions()

    # adjust the selection limits if needed
    if not tryChandra:
        OC.selectionLims[colMag] = [magLo, magHi]
        #OC.selectionLims['FLAGS'] = [0,3]
    OC.selectByLimits()
    
    OC.findAlignedBBox()
    OC.showSelectedPoints()

    OC.getAstCat()
    OC.writeAstCat()

    # Try putting both onto the tangent plane
    # OC.deprojectBothToTangentPlane()
    
    # OC.queryExternal()

    print OC.tAst.colnames
    
    #return

    # try projecting the coords onto the tangent plane
    # READ IN HEADER AND UPDATE
    
    # get CRVALs from fits header? Default to the median RA, DEC from the field

    print "Starting matching..."
    OC.matchOnSphere()

    OC.showDistributions(doLog=showLog)
    return

    OC.showSelectedPoints()

    # nudge the image header
    OC.updateImageRefCoo(img2Shift)

    # write the astrom catalog to disk - this time using our
    # class-level variable
    OC.writeAstCat()
    # OC.writeAstCat('TEST_astromCat.csv')
    
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
