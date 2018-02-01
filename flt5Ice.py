# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:30:28 2017

@author: calorim
"""

import numpy as np
import matplotlib.pyplot as mpl
from scipy import convolve
from scipy.optimize import curve_fit

flt6dat = np.transpose(np.loadtxt('k8r61_spectrum.dat'))
flt6x = (flt6dat[0]+flt6dat[1])/2.
flt6y = flt6dat[2]
flt6eff = 251.*flt6dat[3]

# assumes same binning 
# also assumes file was made without including ice
flt5dat = np.transpose(np.loadtxt('j5r62_spectrum.dat'))
flt5x = (flt5dat[0]+flt5dat[1])/2.
flt5y = flt5dat[2]
flt5eff = 281.*flt5dat[3]

# sum nb bins together
nb = 25
rb6y = np.diff(np.cumsum(flt6y/flt6eff)[::nb])
rb6e = np.sqrt(np.diff(np.cumsum(flt6y/np.square(flt6eff))[::nb]))/rb6y
rb5y = np.diff(np.cumsum(flt5y/flt5eff)[::nb])
rb5e = np.sqrt(np.diff(np.cumsum(flt5y/np.square(flt5eff))[::nb]))/rb5y

# average energy of superbin
rbx = np.diff(np.cumsum(flt6x*(flt6y/flt6eff+flt5y/flt5eff))[::nb])/(rb6y+rb5y)

# throughput-corrected ratio and error
ratio = rb5y/rb6y
error = np.sqrt(np.square(rb6e) + np.square(rb5e))*ratio

gp = (rbx>.2)&(rbx<2)
fity = ratio[gp]
fite = error[gp]
fitx = 1000*rbx[gp]

def iceTrans(e,t,f):
    iceTable = np.loadtxt('data/length/O.len')
    iceRo = 1.
    # convert column depth to thickness in um
    iceAt = (10**4)*np.interp(e,iceTable[:,0],iceTable[:,1])/iceRo
    absorption = 1.-np.exp(-1.*t/iceAt)
    transmission = 1.-f*absorption
    return transmission
    
def avgIce(e,t,f):
    it = iceTrans(1000*flt5x,t,f)
    ait = np.diff(np.cumsum(it*flt5y/flt5eff)[::nb])/rb5y
    return ait[gp] 
    
f = lambda x, *p: avgIce(x,p[0],p[1])
guess = np.ones(2)

popt,pcov = curve_fit(f,fitx,fity,p0=guess,sigma=fite,absolute_sigma=True)

mpl.errorbar(fitx,fity,fite,fmt='.')
mpl.plot(fitx,f(fitx,*popt),
         label='%.1f um Ice, %.0f%% Coverage' % (popt[0],100*popt[1]))
mpl.errorbar(fitx,fity-f(fitx,*popt),fite,fmt='.')
mpl.axhline(0)
mpl.xlabel('Energy (eV)')
mpl.ylabel('Ratio of Smoothed Scaled Counts')
mpl.title('j5r62/k8r61 Contamination')
mpl.legend(loc='best')
mpl.grid()
mpl.show(block=True)

np.savetxt('data/j5r62_ice.dat', popt)
    



