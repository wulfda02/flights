# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:13:28 2017
@author: dwulf
"""

import flightLibrary as fl
import warnings
import os

warnings.filterwarnings("ignore",".*GUI is implemented.*")
warnings.filterwarnings("ignore",".*Covariance of the parameters.*")
warnings.filterwarnings("ignore",".*invalid value encountered*")
warnings.filterwarnings("ignore",".*Polyfit may be poorly conditioned")

run = 'k8r61'
pixels = range(36)

# Make bsn file from card data if not already done
# Also get timing
fn = "%s.bsn.hdf5" % run
if not os.path.isfile(fn):
    
    print "Making bsn file..."
    cardData = fl.cardData(run)
    cardData.writeAllPixels()
    
    print "Applying timing..."
    timing = fl.groupTiming(run)
    timing.correctedFlightTiming()
    

#for pixel in pixels:
#    print run, pixel
#
#    # Load data for pixel from hdf5 file
#    pxlObj = fl.singlePixel(run,pixel)
#    pxlObj.loadFromHdf('bsn')
#    
#    # Make optimum filters
#    filters = pxlObj.getFilter()
#    if filters is None:
#        print "Bad Pixel: Unable To Make Filter"
#        continue
#    if (filters.bslnRes()>15)and(pixel!=18):
#        print "Bad Pixel: Poor Resolution"
#        continue
#        
#    # Convolve filters with data
#    pxlObj.filterData(filters)
#    
#    # Fit filtered data
#    pxlObj.fitFilteredData(filters)
#    
#    # Save fit to hdf5 file
#    pxlObj.writeHdfFit('bsn')
#        
## Apply temperature-corrected non-linear gain scale
runObj = fl.allPixels(run)
#runObj.eScale()

# Create on-target spectrum
obsObj = runObj.selectEvents()
phaFile = obsObj.spectrum([50,2000],2.5)

import xspec as xs
xs.Spectrum(phaFile)
#xs.Model("vmekal+wabs*(vmekal+bknpower+bknpower)")
#xs.Fit.perform()
#xs.Plot.device = "/xs"
#xs.Plot("data")



