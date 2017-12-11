# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 13:46:34 2017

@author: calorim
"""

import xspec as xs
import numpy as np
import matplotlib.pyplot as mpl

# load spectrum
s = xs.Spectrum("k8r61.pha")

# solar abundance (Anders and Grevesse 1989)
# cool cloud depleted abundances (Savage and Sembach 1996)
xs.Xset.abund = "angr"
#xs.Xset.abund = "file sscooln.abund"

cmPerPc = 3.086e+18 # pc to cm conversion

## unabsorbed thermal plus absorbed thermal and powerlaw
#m = xs.Model("apec + wabs*(apec + pow)")
## set initial parameters to McCammon 2002 values (solar)
#m.setPars(.099,1.,0.,.0088*1e-14*cmPerPc/(4*np.pi),
#          .018,
#          .225,1.,0.,.0037*1e-14*cmPerPc/(4*np.pi),
#          1.52,12.3)
## freeze power law
#m.powerlaw.PhoIndex.frozen = True
#m.powerlaw.norm.frozen = True
          
# unabsorbed thermal plus absorbed thermal plus
# broken powerlawfrom Gwynne's thesis
# set initial plasma to McCammon 2002 values (solar)
m = xs.Model("apec + wabs*(apec + bknp + bknp)")

m.setPars(.099,1.,0.,.0088*1e-14*cmPerPc/(4*np.pi),
          .018,
          .225,1.,0.,.0037*1e-14*cmPerPc/(4*np.pi),
          1.54,.9,1.4,5.7,
          1.96,.9,1.4,4.9)
# freeze power laws
m.bknpower.PhoIndx1.frozen = True
m.bknpower.BreakE.frozen = True
m.bknpower.PhoIndx1.frozen = True
m.bknpower_5.norm.frozen = True
m.bknpower_5.PhoIndx1.frozen = True
m.bknpower_5.BreakE.frozen = True
m.bknpower_5.PhoIndx1.frozen = True
m.bknpower_5.norm.frozen = True

# freeze absorbtion
m.wabs.nH.frozen = True

# scheme for estimating variance (needed for empty bins)
xs.Fit.weight = "model"
#xs.Fit.statMethod = "cstat"
s.ignore("**-.2 1.-**")

xs.Fit.perform()

xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")

energy = 1000*np.array(xs.Plot.x())
counts = np.array(xs.Plot.y())
model = np.array(xs.Plot.model())
mpl.step(energy,counts,where="mid",lw=.5)
mpl.plot(energy,model)
mpl.xlabel("Energy (eV)")
mpl.ylabel("cts/sec/bin")
mpl.grid(ls=":")

def unityGauss(x,e,fwhm):
    sig = fwhm/2.355
    exp = -0.5*np.square((x-e)/sig)
    return np.exp(exp)/(np.sqrt(2*np.pi)*sig)
    
def unityFlat(x,e):
    y = np.zeros(len(x))
    i = np.where(x<e)
    y[i] = 1./e
    return y
    
def peRsp(x,e,fwhm,f):
    return f*unityFlat(x,e) + (1.-f)*unityGauss(x,e,fwhm)
    
x = np.arange(2000)
y = peRsp(x,1500,10.,.05)
mpl.plot(x,y)
mpl.show(block=True)
