# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:13:28 2017
@author: dwulf
"""

import flightLibrary as fl
import warnings
import sys
import os

warnings.filterwarnings("ignore",".*GUI is implemented.*")
warnings.filterwarnings("ignore",".*Covariance of the parameters.*")
warnings.filterwarnings("ignore",".*invalid value encountered*")

# Uses same options/logic as fitsfilter
if len(sys.argv)>1:
    run = sys.argv[1]
    if len(sys.argv)>2:
        pixel = int(sys.argv[2])
    else:
        pixel = -1
        
    if len(sys.argv)>3 and pixel>=0:
        autoFlag = int(sys.argv[3])
    elif pixel>=0:
        autoFlag = 0
    else:
        autoFlag = 1
else:
    scriptName = os.path.basename(__file__)
    print "Usage: python %s runBaseName [pixel [autoflag]]" % scriptName
    print "runBaseName e.g. k8r61"
    print "pixel 0-35 (-1 for batch mode [default])"
    print "autoflag 0 for single pixel interactive"
    print "autoFlag 1 for otherwise"
    run = 'k8r61'
    pixel = 0
    autoFlag = 0
    #sys.exit()
    
if pixel<0:
    pixels = range(35)
else:
    pixels = [pixel]
interactive = not autoFlag

# Make bsn file from card data if not already done
# Also get timing
fn = "%s.bsn.hdf5" % run
if not os.path.isfile(fn):
    print "Making bsn file..."
    cardData = fl.cardData(run)
    cardData.writeAllPixels()
    print "Getting timing from bias toggles..."
    timing = fl.groupTiming(run)
    timing.correctedFlightTiming()
    

for pixel in pixels:
    print run, pixel

    # Load data for pixel from hdf5 file
    pxlObj = fl.singlePixel(run,pixel,interactive)
    pxlObj.loadFromHdf('bsn')
    
    # Make optimum filters
    filters = pxlObj.getFilter()
    if filters is None:
        print "Bad Pixel: Unable To Make Filter"
        continue
    if (filters.bslnRes()>15) and not interactive:
        print "Bad Pixel: Poor Resolution"
        continue
        
    # Convolve filters with data
    pxlObj.filterData(filters)
    
    # Fit filtered data
    pxlObj.fitFilteredData(filters)
    
    # Apply non-linear gain scale
    pxlObj.nonLinGain()
    
    # Save fit to hdf5 file
    save = True
    if interactive:
        save = raw_input("Save session? (Y/n)")[:1].lower()=='y'
    if save:
        pxlObj.writeHdfFit('bsn')

