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

## Get timing from swp file    
#print "Applying timing..."
#timing = fl.groupTiming(run)
#timing.swpTiming()
#    
#
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
        
# Apply temperature-corrected non-linear gain scale
runObj = fl.allPixels(run)
#runObj.eScale()

import matplotlib.pyplot as mpl
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

def unityFlat(x,e):
    y = np.zeros(x.shape)
    i = np.where(x<e)
    y[i] = 1./e
    return y
    
def unityGauss(x,e,fwhm):
    sig = fwhm/2.355
    exp = -0.5*np.square((x-e)/sig)
    return np.exp(exp)/(np.sqrt(2*np.pi)*sig)
    
def unityGaussWTail(x,e,fwhm,tau): # emperical microcalorimeter line shape
    sig = fwhm/2.355   
    return (1./(2.*tau)
            *np.exp((x-e)/tau + sig**2/(2.*tau**2))
            *erfc((1/np.sqrt(2.))*((x-e)/sig + sig/tau)))

def unityLorentz(x,e,fwhm):
    return ((1./np.pi)*(0.5*fwhm/(np.square(x-e)+np.square(0.5*fwhm))))
    
def unityTwoGauss(x,e,fw1,fw2,de,f):
    return (1.-f)*unityGauss(x,e,fw1) + f*unityGauss(x,e+de,fw2)
    
def unityGaussLorentz(x,e,fw1,fw2,f):
    return (1.-f)*unityGauss(x,e,fw1) + f*unityLorentz(x,e,fw2)

    
absDat = np.transpose(np.loadtxt('data/absorberTransmission.dat'))
siDat = np.transpose(np.loadtxt('data/substrateTransmission.dat'))    

guess = [3313,200000,
         3589,20000,
         15,15,.1,
         .83,150,
         .05,1,
         3]

def fitCal(x,KaE,KaA,KbE,KbA,FW1,FW2,dF,MdE,MFWHM,peF,cor,bg):
    KaAbsT = np.interp(KaE,absDat[0],absDat[1])+.03
    KaSubT = np.interp(KaE,siDat[0],siDat[1])
    KaMushF = cor*(KaAbsT*(1.-KaSubT))/(1.-KaAbsT)
    Ka = (1.-peF)*KaA*((1.-KaMushF)*unityGaussLorentz(x,KaE,FW1,FW2,dF)
                       +KaMushF*unityGauss(x,MdE*KaE,MFWHM))
    KbAbsT = np.interp(KbE,absDat[0],absDat[1])+.03
    KbSubT = np.interp(KbE,siDat[0],siDat[1])
    KbMushF = cor*(KbAbsT*(1.-KbSubT))/(1.-KbAbsT)
    Kb = (1.-peF)*KbA*((1.-KbMushF)*unityGaussLorentz(x,KbE,FW1,FW2,dF)
                       +KbMushF*unityGauss(x,MdE*KbE,MFWHM*KbE/KaE))
    KaPe = peF*KaA*unityFlat(x,KaE)
    KbPe = peF*KbA*unityFlat(x,KbE)
    return Ka + Kb + KaPe + KbPe + bg

pls = runObj.pulses()
plsEnrg = pls['energy']
plsTm = pls['time']
plsPxl = pls['pixel']
#mpl.scatter(plsTm,plsEnrg,marker="+",lw=.2)
#mpl.show(block=True)

lastGoodTime = 88500.
badPixels = [2,7,8,9,15,18,19,21,24,26,27,29,31,35]
rng = [1800,5000]
bnw = 10.
bns = int((rng[1]-rng[0])/bnw)
Hy = np.zeros(bns)
for p in pixels:
    if p not in badPixels:
        print p
        goodPls = np.where((plsTm<lastGoodTime)&(plsPxl==p))
        pltPls = plsEnrg[goodPls]       
        hy,hx =  np.histogram(pltPls,bins=bns,range=rng)
        Hx = hx[:-1] + 0.5*bnw
        Hy += hy
        
He = np.sqrt(Hy)  
He[np.where(He==0)] = 1

f = lambda x,*p: bnw*fitCal(x,*p)

popt,pcov = curve_fit(f,Hx,Hy,p0=guess,sigma=He,absolute_sigma=True,bounds=(0,np.inf))
perr = np.sqrt(np.diag(pcov))
        
r = Hy-f(Hx,*popt)
mpl.figure()
mpl.subplot(2,1,1)
mpl.step(Hx,Hy,where='mid')
mpl.step(Hx,f(Hx,*popt),where='mid')
mpl.subplot(2,1,2)
mpl.errorbar(Hx,r,He,fmt=".")
mpl.axhline(0)
mpl.show(block=True)






