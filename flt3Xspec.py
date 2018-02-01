# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 13:46:34 2017

@author: calorim
"""

import xspec as xs
import numpy as np
import matplotlib.pyplot as mpl

# define own absorption model
# based on IR nH data, using 
# Ander & Ebihara metalicities
# Wisconsin cross sections
# and averaged over FOV
fn = "data/irnh_AndersEbihara_90-60.dat"
absorption = np.transpose(np.loadtxt(fn))
keV = .001*absorption[0]
nH = absorption[1]
sigma = absorption[2]

def avgabs(energies,params,factor):
    for i in range(len(energies)-1):
        eLo = energies[i]
        eHi = energies[i+1]
        nSteps = max(2,int(1000*(eHi-eLo)+1)) # approx 1 eV steps 
        e = np.linspace(eLo,eHi,nSteps)
        col = params[0]*np.interp(e,keV,nH)
        sig = np.interp(e,keV,sigma)
        avgT = np.trapz(np.exp(-1*col*sig),e)/(eHi-eLo)
        factor[i] = avgT

absInfo = ("norm '' 1.0 0.0 0.1 1.9 2.0 0.00001",)
xs.AllModels.addPyMod(avgabs,absInfo,'mul',spectrumDependent=True)  

# Add Randall Smith's CX model
xs.AllModels.lmod("acx","/home/calorim/Downloads/acx-1.0.1/xspec")    

# load spectrum
s = xs.Spectrum("flt3.pha")

# solar abundance (Anders and Grevesse 1989)
xs.Xset.abund = "angr"
z = "Solar"
# cool cloud depleted abundances (Savage and Sembach 1996)
#xs.Xset.abund = "file sscooln.abund"
#z = "Depleted"

cmPerPc = 3.086e+18 # pc to cm conversion

# absorbed broken powerlawfrom Mushotzky 2000 (AGN)
# + absorbed thermal (Halo)
# + unabsorbed thermal (LHB)
# + cal lines
# thermal norms initially set to 0
m = xs.Model("apec + avgabs*(bknp + bknp + apec)")
m.setPars(.08617,1,0,0,
          1.,
          1.54,1.2,1.4,5.7,
          1.96,1.2,1.4,4.9,
          .2,1,0,0)
          
# freeze power laws
m.bknpower.PhoIndx1.frozen = True
m.bknpower.BreakE.frozen = True
m.bknpower.PhoIndx2.frozen = True
m.bknpower.norm.frozen = True
m.bknpower_4.PhoIndx1.frozen = True
m.bknpower_4.BreakE.frozen = True
m.bknpower_4.PhoIndx2.frozen = True
m.bknpower_4.norm.frozen = True

# freeze absorbtion
m.avgabs.norm.frozen = True

# freeze LHB T
#m.apec.kT.frozen = True

# scheme for estimating variance (needed for empty bins)
xs.Fit.weight = "model"

# set fit statistic
xs.Fit.statMethod = "pgstat"

xs.Fit.perform()

# access important numbers
kb = 8.617e-8 # keV/K
LHBT = m.apec.kT.values[0]/kb
LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
HaloT = m.apec_5.kT.values[0]/kb
HaloEM = m.apec_5.norm.values[0]*4*np.pi/(1e-14*cmPerPc)

# set plotting parameters
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
mpl.text(200,20,"LHB T: %.1e K\nLHB EM: %.1e pc cm^-6" % (LHBT,LHBEM))
mpl.text(700,30,"Halo T: %.1e K\nHalo EM: %.1e pc cm^-6" % (HaloT,HaloEM))
mpl.title("%s Abundances" % z)
mpl.grid(ls=":")
mpl.show(block=True)




