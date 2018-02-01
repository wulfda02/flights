# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:30:28 2017

@author: calorim
"""

import numpy as np
import matplotlib.pyplot as mpl
import xspec as xs

flt6dat = np.transpose(np.loadtxt('k8r61_spectrum.dat'))
flt6datx = (flt6dat[0]+flt6dat[1])/2.
flt6daty = flt6dat[2]
flt6eff = 251.*flt6dat[3]
flt6src = flt6daty/flt6eff
flt6err = np.sqrt(flt6daty)
flt6err[np.where(flt6err==0)] = 1.
flt6err = flt6err/flt6eff

s = xs.Spectrum("k8r61.pha")
xs.Xset.abund = "wilm"
xs.Xset.xsect = "vern"
sr = 0.809259 # XQC FOV in sr
m = xs.Model("phabs*(bknp)") # Kirsch
m.setPars(.407,2.02,2.,2.07,8.95*sr)
xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")
crabobs = np.array(xs.Plot.model())
crabsrc = crabobs/flt6eff

rosatresp = np.loadtxt('data/rosat_rsp.dat',skiprows=1,usecols=(0,1,2,3,4,5,6,7)) # ROSAT R124567 response
binCenter,binWidth,r1,r2,r4,r5,r6,r7= np.transpose(rosatresp)
binCenter = binCenter
binWidth = binWidth

bandDict = {1:r1,2:r2,4:r4,5:r5,6:r6,7:r7}
rosatCts = {1:192.46,2:205.96,4:32.91,5:49.65,6:70.16,7:43.71}

mpl.step(flt6datx,flt6daty,where='mid')
mpl.plot(flt6datx,crabobs)
for band in bandDict:
    rosatSpec6 = flt6src*np.interp(flt6datx,binCenter,bandDict[band])
    rosatSpec6Err = flt6err*np.interp(flt6datx,binCenter,bandDict[band])
    rosatCrab = crabsrc*np.interp(flt6datx,binCenter,bandDict[band])
    brightness6 = np.sum(rosatSpec6)
    brightness6Err = np.sqrt(np.sum(np.square(rosatSpec6Err)))
    brightnessCrab = np.sum(rosatCrab)
    netBright = brightness6 - brightnessCrab
    print "%d: %.1f +/- %.1f (%.1f)" % (band,netBright,brightness6Err,rosatCts[band])
    mpl.plot(binCenter,bandDict[band])
mpl.show(block=True)





