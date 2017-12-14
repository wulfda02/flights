# -*- coding: utf-8 -*-
"""
Created on Dec 14 2017
Used to determine calibration source
photo-electric escape background.
This is not a true empty belljar, but
it should be ok.
@author: dwulf

"""

import flightLibrary as fl
import warnings
import os

warnings.filterwarnings("ignore",".*GUI is implemented.*")
warnings.filterwarnings("ignore",".*Covariance of the parameters.*")
warnings.filterwarnings("ignore",".*invalid value encountered*")
warnings.filterwarnings("ignore",".*Polyfit may be poorly conditioned")
warnings.filterwarnings("ignore",".*divide by zero encountered*")

run = 'l5r18'
pixels = range(36)

# Make bsn file from card data if not already done
fn = "%s.bsn.hdf5" % run
if not os.path.isfile(fn):
    
    print "Making bsn file..."
    cardData = fl.cardData(run)
    cardData.writeAllPixels()

# Get timing from swp file    
print "Applying timing..."
timing = fl.groupTiming(run)
timing.swpTiming()
    

for pixel in pixels:
    print run, pixel

    # Load data for pixel from hdf5 file
    pxlObj = fl.singlePixel(run,pixel)
    pxlObj.loadFromHdf('bsn')
    
    # Make optimum filters
    filters = pxlObj.getFilter()
    if filters is None:
        print "Bad Pixel: Unable To Make Filter"
        continue
    if (filters.bslnRes()>15)and(pixel!=18):
        print "Bad Pixel: Poor Resolution"
        continue
        
    # Convolve filters with data
    pxlObj.filterData(filters)
    
    # Fit filtered data
    pxlObj.fitFilteredData(filters)
    
    # Save fit to hdf5 file
    pxlObj.writeHdfFit('bsn')
        
## Apply temperature-corrected non-linear gain scale
#runObj = fl.allPixels(run)
#runObj.eScale()

## IR filter stack
#filters = {'IDB':['11895-2.dat','polyimide','none'],
#           'ODB':['ODB103_092013.dat','polyimide','none'],
#           'I2K':['2K_30048-3.dat','polyimide','new'],
#           'O2K':['2K105_082013.dat','polyimide','new'],
#           '130K':['130K116_082013.dat','polyimide','new'],
#           'RT':['RT114_092013.dat','polyimide','new']}
#filterStack = open("data/%s_filterstack.dat" % run,"w")
#for fid in filters:
#    f = filters[fid]
#    irf = fl.IRFilter()
#    filterStack.write(irf.fitFromFile(*f))
#filterStack.close()

# Create on-target spectrum
#obsObj = runObj.selectEvents()
# specra with significantly more than 1000 bins may cause memory 
# problems when making the response matrix
#phaFile = obsObj.spectrum([100,4000],3,fileNameAdd="test")

# Xspec fitting
#import xspec as xs
#xs.Spectrum(phaFile)
#xs.Model("vmekal+wabs*(vmekal+bknpower+bknpower)")
#xs.Fit.perform()
#xs.Plot.device = "/xs"
#xs.Plot("data")



