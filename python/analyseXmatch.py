#
# analyseXmatch.py
#

# WIC & VF 2017-05-19 - some methods to analyse the cross-matches to
# Annalisa Calamida's stroemgren photomertry

########## IMPORT STATEMENTS ############

# Useful to read in the data
import os
from astropy.table import Table

# We'll almost certainly need this to do anything interesting:
import numpy as np

# do we have fancy plotting?
try:
    import seaborn as sns
    HAVESEABORN=True
except:
    HAVESEABORN=False

# for cutting down our sample
from astropy.stats import sigma_clip

# it's useful to plot what we have...
import matplotlib.pylab as plt

# configure pylab to allow returning to the interpreter after plotting.
plt.ion()  
plt.style.use('ggplot') # make our plots look fancy

################ Import stuff finished. Now for the methods. #####

def loadTable(pathTab='', Verbose=True):

    """Loads the specified table to memory, returns an astropy table object"""

    # Ensure we always return a table, even if it is blank
    tRet = Table()

    # some defensive programming: if no path specified...
    if len(pathTab) < 2:
        if Verbose:
            print "analyseXmatch.loadTable FATAL - input path not specified." 
        return tRet

    if not os.access(pathTab, os.R_OK):
        if Verbose:
            print "analyseXmatch.loadTable FATAL - cannot read path %s" \
                % (pathTab)
        return tRet

    try:
        tRet = Table.read(pathTab, format='ascii')
    except:
        print "analyseXmatch.loadTable WARN - problem reading %s" \
            % (pathTab)

    # whatever's happened, by this point we have at least a blank
    # table. So return it here.
    return tRet

def selectOnMagnitude(tPhot=Table(), magMax=22., deltaMax=4.0, \
                          colMag='Instrumental_VEGAMAG', \
                          colDel='A_Magn'):

    """Returns boolean selection indices for objects satisfying our
    selection region. Currently the selection boxes are hard-coded.""" 

    # (we could update later to do something clever when selecting the
    # regions automatically...)
    
    if np.sum(checkColumns(tPhot, [colMag, colDel])) < 2:
        print "selectOnMagnitude WARN - at least one column missing."
        return

    mag = tPhot[colMag]
    delta = tPhot[colDel]

    # let's do our boolean now
    bGood = (mag <= magMax) & (delta <= deltaMax)

    return np.asarray(bGood)

def calcStats(data=[], nIters=4, nSigma=3):    

    """Calculate some straightforward statistics for the subsample of
    interest."""

    # This is where more sophisticated methods might come in...

    # initialize empty "dictionary" of statistics.
    DStats = {}
    
    if len(data) < 1:
        return DStats

    # Let's just write out the stats we want for the moment... we can
    # refine this later.
    DStats['median'] = np.median(data)
    DStats['stddev'] = np.std(data)
    DStats['nPoints'] = np.size(data)

    # let's try something more interesting...  array of objects that
    # survived sigma_clipping
    dataClipped = sigma_clip(data, sigma=nSigma, iters=nIters)

    DStats['sigClip_mean'] = np.mean(dataClipped)
    DStats['sigClip_std'] = np.std(dataClipped)
    DStats['sigClip_median'] = np.median(dataClipped)

    # doesn't hurt to also send back up the clipping parameters
    DStats['sigClip_params'] = {'nIters': nIters, 'nSigma':nSigma}
    DStats['sigClip_params']['nSurvived'] = np.sum(~dataClipped.mask)

    return DStats

def showDeltaMags(tPhot=Table(), \
                      colMag='Instrumental_VEGAMAG', \
                      colDel='A_Magn', \
                      figNum=1, \
                      bSel=np.array([]), \
                      showFilename=True, \
                      figName='./TEST_deltaMags.png'):

    """Plots the delta-magnitudes. Overplots "good" magnitudes in a
    different color."""

    # do our columns exist?
    if not colMag in tPhot.colnames:
        print "analyseXmatch.showDeltaMags WARN - mag col %s not in table." \
            % (colMag)
        return

    if not colDel in tPhot.colnames:
        print "analyseXmatch.showDeltaMags WARN - delta col %s not in table." \
            % (colDel)
        return

    # if we are here, then we have our deltas and can
    # proceed. Convenience views:
    mag = tPhot[colMag]
    delta = tPhot[colDel]

    fig = plt.figure(figNum)
    fig.clf()

    ax1=fig.add_subplot(111)
    ax1.scatter(mag, delta, s=5, alpha=0.5)
    
    ax1.set_xlabel(colMag)
    ax1.set_ylabel(colDel)

    # if the table has the 'filename' attribute, show it.
    if 'filename' in tPhot.meta.keys() and showFilename:
        sName = tPhot.meta['filename']
        ax1.annotate(sName, (0.05,0.97), xycoords='axes fraction', \
                         ha='left', va='top')

    # that's the whole lot... let's overplot the "good" objects too.
    if np.sum(bSel) > 0:
        ax1.scatter(mag[bSel], delta[bSel], color='g')
    
    fig.savefig(figName)

def showHistSubset(tSubset=Table(), colDel='A_Magn', figNum=2, \
                       nBins=100, xRange=(0,3), figName='./TEST_histo.png'):

    """Quick routine to show the histogram of the input table"""

    if not colDel in tSubset.colnames:
        return

    if len(tSubset) < 1:
        return

    fig = plt.figure(figNum)
    fig.clf()

    ax1=fig.add_subplot(111)
    
    # data to plot
    data = tSubset[colDel]

    dum = ax1.hist(data, bins=nBins, alpha=0.75, range=xRange)

    ax1.set_xlabel(colDel)

    fig.savefig(figName)

def checkColumns(tPhot=Table(), LCols=[] ):

    """Tests if all of the columns in list LCols are present in table
    tPhot."""

    LFound = [] # initialise

    if len(LCols) < 1:
        return LFound

    # Go through each column in turn, test if present.
    for iCol in range(len(LCols)):
        if LCols[iCol] in tPhot.colnames:
            LFound.append(True)
        else:
            LFound.append(False)

    return LFound

######################### Test routines follow ##################

def testReadTable(pathMatch='./F410M-Calamida.dat', \
                      magMax=23., delMax=3.5, \
                      colDelt='A_Magn', Verbose=True):

    """Tests our routines. pathMatch = path to the matched
    table. Examples (on my system I renamed F467M-Calamida to
    F467M-Calamida.dat, but it's probably not necessary):

    analyseXmatch.testReadTable()

    analyseXmatch.testReadTable('./F467M-Calamida.dat', colDelt='My_data_-_Calamida', magMax=21.6, Verbose=False)"""

    # WIC on my laptop, the files happen to lie in 
    # /Users/clarkson/Data/stromgren/matched2AC

    tPhot = loadTable(pathMatch)

    # just for convenience, let's add the filename to the table
    # metadata so that we can access it later:
    tPhot.meta['filename'] = pathMatch[:]
    
    # let's see what we have.
    if Verbose:
        print tPhot.colnames

    # test our column checker
    bAccept = selectOnMagnitude(tPhot, magMax, delMax, \
                                    colDel=colDelt)

    # let's compute some statistics!
    DChar = calcStats(tPhot[colDelt][bAccept])

    # just dump to screen for the moment
    for sKey in DChar.keys():
        print sKey, DChar[sKey]

    # let's take a look at the deltas... what do we have here?
    showDeltaMags(tPhot, colDel=colDelt, bSel=bAccept)

    # Also show the histogram, just for information (there are other
    # ways to visualize such small samples...) Consider seaborn.rugplot...
    showHistSubset(tPhot[bAccept], colDelt)

    # let's build a rugplot here: (not sure why this isn't working...)
    if HAVESEABORN:
        plt.figure(3, figsize=(8,5))
        plt.clf()

        #g = np.where(bAccept)
        #vec = np.asarray(tPhot[g])

        sns.distplot(tPhot[bAccept][colDelt], hist=True, rug=True, color='0.2')
        plt.title('KDE of subset of deltas')
        plt.ylabel('KDE (%s)' % (colDelt))
        plt.xlabel(colDelt)

        plt.savefig('./TEST_rugplot.png')
