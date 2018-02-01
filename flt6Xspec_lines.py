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
# Ander & Ebihara abundances
# Wisconsin cross sections
# and averaged over FOV
fn = "data/irnh_AndersEbihara_165--5.dat"
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

# load spectrum
s = xs.Spectrum("k8r61.pha")

xs.Xset.abund = "angr"
xs.Xset.xsect = "vern"

cmPerPc = 3.086e+18 # pc to cm conversion
sr = 0.809259 # XQC FOV in sr

agn = "avgabs*(bknp + bknp)" # Mushotzky 2000
crab = "vphabs*(bknp)" # Kirsch using wilm abundances
cal = "gau + gau"
halo = "avgabs*(vapec)" # Savage sembach abund, no C,O,Ne
lhb = "vapec" # Savage sembach abund, no C,O,Ne
ovii = "gau + gau + gau" # Ka triplet
oviii = "gau + gau + gau" # Ly-a,b,c
neix = "gau + gau + gau" # Ka triplet
nex = "gau + gau + gau" # Ly-a,b,c
cvi = "gau + gau + gau" # Ly-a,b,c
cv = "gau + gau + gau + gau" # Ka triplet + Kb
mod = " + ".join([agn, crab, cal, halo, lhb, ovii, oviii, neix, nex, cvi, cv])
m = xs.Model(mod)
m.setPars(1.,
          1.54,1.2,1.4,5.7,
          1.96,1.2,1.4,4.9,
          .407,.81,.67,.60,.58,.60,.93,.93,.98,.98,.96,.58,.98,.96,.96,.63,.98,.93,
          2.02,2.,2.07,8.95*sr,
          3.318,.01,25.,
          3.589,.01,5.,
          1.,
          .17,1,0,.71,0,0,.028,1,.049,1.74,.33,.00018,.0037,.93,0,0,
          .085,1,0,.71,0,0,.028,1,.049,1.74,.33,.00018,.0037,.93,0,0,
          .5613,0,0,
          .5689,0,0,
          .5748,0,0,
          .6536,0,0,
          .7746,0,0,
          .8171,0,0,
          .9055,0,0,
          .9155,0,0,
          .9231,0,0,
          1.0217,0,0,
          1.2112,0,0,
          1.2773,0,0,
          .3675,0,0,
          .4355,0,0,
          .4594,0,0,
          .299,0,0,
          .3045,0,0,
          .3085,0,0,
          .3538,0,0)
          
# freeze power laws
m.bknpower.PhoIndx1.frozen = True
m.bknpower.BreakE.frozen = True
m.bknpower.PhoIndx2.frozen = True
m.bknpower.norm.frozen = True
m.bknpower_3.PhoIndx1.frozen = True
m.bknpower_3.BreakE.frozen = True
m.bknpower_3.PhoIndx2.frozen = True
m.bknpower_3.norm.frozen = True
m.bknpower_5.PhoIndx1.frozen = True
m.bknpower_5.BreakE.frozen = True
m.bknpower_5.PhoIndx2.frozen = True
m.bknpower_5.norm.frozen = True

# freeze absorbtion
m.avgabs.norm.frozen = True
m.vphabs.nH.frozen = True
m.avgabs_8.norm.frozen = True

# freeze halo
m.vapec.norm.frozen = True
m.vapec.kT.frozen = True

# freeze LHB
m.vapec_10.norm.frozen = True
m.vapec_10.kT.frozen = True

# freeze ion lines
m.gaussian_11.LineE.frozen = True
m.gaussian_11.Sigma.frozen = True
m.gaussian_11.norm.frozen = True
m.gaussian_12.LineE.frozen = True
m.gaussian_12.Sigma.frozen = True
m.gaussian_12.norm.frozen = True
m.gaussian_13.LineE.frozen = True
m.gaussian_13.Sigma.frozen = True
m.gaussian_13.norm.frozen = True
m.gaussian_14.LineE.frozen = True
m.gaussian_14.Sigma.frozen = True
m.gaussian_14.norm.frozen = True
m.gaussian_15.LineE.frozen = True
m.gaussian_15.Sigma.frozen = True
m.gaussian_15.norm.frozen = True
m.gaussian_16.LineE.frozen = True
m.gaussian_16.Sigma.frozen = True
m.gaussian_16.norm.frozen = True
m.gaussian_17.LineE.frozen = True
m.gaussian_17.Sigma.frozen = True
m.gaussian_17.norm.frozen = True
m.gaussian_18.LineE.frozen = True
m.gaussian_18.Sigma.frozen = True
m.gaussian_18.norm.frozen = True
m.gaussian_19.LineE.frozen = True
m.gaussian_19.Sigma.frozen = True
m.gaussian_19.norm.frozen = True
m.gaussian_20.LineE.frozen = True
m.gaussian_20.Sigma.frozen = True
m.gaussian_20.norm.frozen = True
m.gaussian_21.LineE.frozen = True
m.gaussian_21.Sigma.frozen = True
m.gaussian_21.norm.frozen = True
m.gaussian_22.LineE.frozen = True
m.gaussian_22.Sigma.frozen = True
m.gaussian_22.norm.frozen = True
m.gaussian_23.LineE.frozen = True
m.gaussian_23.Sigma.frozen = True
m.gaussian_23.norm.frozen = True
m.gaussian_24.LineE.frozen = True
m.gaussian_24.Sigma.frozen = True
m.gaussian_24.norm.frozen = True
m.gaussian_25.LineE.frozen = True
m.gaussian_25.Sigma.frozen = True
m.gaussian_25.norm.frozen = True
m.gaussian_26.LineE.frozen = True
m.gaussian_26.Sigma.frozen = True
m.gaussian_26.norm.frozen = True
m.gaussian_27.LineE.frozen = True
m.gaussian_27.Sigma.frozen = True
m.gaussian_27.norm.frozen = True
m.gaussian_28.LineE.frozen = True
m.gaussian_28.Sigma.frozen = True
m.gaussian_28.norm.frozen = True
m.gaussian_29.LineE.frozen = True
m.gaussian_29.Sigma.frozen = True
m.gaussian_29.norm.frozen = True

# set gaussian widths equal
m.gaussian_7.Sigma.link = "33"

# scheme for estimating variance (needed for empty bins)
xs.Fit.weight = "model"

# set fit statistic
xs.Fit.statMethod = "cstat"

# ignore low energy bins
s.ignore("**-3.2")

xs.Fit.perform()

# freeze cal source
m.gaussian.LineE.frozen = True
m.gaussian.Sigma.frozen = True
m.gaussian.norm.frozen = True
m.gaussian_7.LineE.frozen = True
m.gaussian_7.Sigma.frozen = True
m.gaussian_7.norm.frozen = True

s.notice("all")
s.ignore("**-.2 1.5-**")

# thaw halo
m.vapec.norm.frozen = False
#m.vapec.kT.frozen = False

# thaw LHB
m.vapec_10.norm.frozen = False

# thaw N abundance
m.vapec.N.frozen = False
m.vapec_10.N.frozen = False

# thaw ion lines
m.gaussian_11.norm.frozen = False
m.gaussian_12.norm.frozen = False
m.gaussian_13.norm.frozen = False
m.gaussian_14.norm.frozen = False
m.gaussian_15.norm.frozen = False
m.gaussian_16.norm.frozen = False
m.gaussian_17.norm.frozen = False
m.gaussian_18.norm.frozen = False
m.gaussian_19.norm.frozen = False
m.gaussian_20.norm.frozen = False
m.gaussian_21.norm.frozen = False
m.gaussian_22.norm.frozen = False
m.gaussian_23.norm.frozen = False
m.gaussian_24.norm.frozen = False
m.gaussian_25.norm.frozen = False
m.gaussian_26.norm.frozen = False
m.gaussian_27.norm.frozen = False
m.gaussian_28.norm.frozen = False
m.gaussian_29.norm.frozen = False

xs.Fit.perform()

# set plotting parameters
xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")

energy = 1000*np.array(xs.Plot.x())
counts = np.array(xs.Plot.y())
model = np.array(xs.Plot.model())

#energy = (energy[:-1:2]+energy[1::2])/2.
#counts = counts[:-1:2]+counts[1::2]
#model = model[:-1:2]+model[1::2]

mpl.fill_between(energy,counts,step='mid',alpha=.4)
mpl.step(energy,counts,where="mid",lw=.5)
mpl.plot(energy,model,lw=2)
mpl.ylabel("cts/sec/bin")
mpl.ylim(ymin=0)
mpl.grid(ls=":")
mpl.xlabel("Energy (eV)")
mpl.show(block=True)




