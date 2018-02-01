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

# Add Randall Smith's CX model
xs.AllModels.lmod("acx","/home/calorim/Downloads/acx-1.0.1/xspec")   

# load spectrum
s = xs.Spectrum("k8r61.pha")

xs.Xset.abund = "angr"
xs.Xset.xsect = "vern"

cmPerPc = 3.086e+18 # pc to cm conversion
sr = 0.809259 # XQC FOV in sr

agn = "avgabs*(bknp + bknp)" # Mushotzky 2000
crab = "vphabs*(bknp)" # Kirsch using wilm abundances
cal = "gau + gau"
halo = "avgabs*(vapec)" # Savage sembach abund
lhb = "vapec" # Savage sembach abund
o = "vacx"
ne = "vacx"
c = "vacx"
fe = "vacx"
mod = " + ".join([agn, crab, cal, halo, lhb, o, ne, c, fe])
m = xs.Model(mod)
m.setPars(1.,
          1.54,1.2,1.4,5.7,
          1.96,1.2,1.4,4.9,
          .407,.81,.67,.60,.58,.60,.93,.93,.98,.98,.96,.58,.98,.96,.96,.63,.98,.93,
          2.02,2.,2.07,8.95*sr,
          3.318,.01,25.,
          3.589,.01,5.,
          1.,
          .17,1,.38,.71,.35,1,.028,1,.049,1.74,.33,.00018,.0037,.93,0,0,
          .085,1,.38,.71,.35,1,.028,1,.049,1.74,.33,.00018,.0037,.93,0,0,
          .17,.1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,8,0,
          .17,.1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,8,0,
          .085,.1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,8,0,
          .17,.1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,8,0)
          
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

# freeze CX
m.vacx.kT.frozen = True
m.vacx.redshift.frozen = True
m.vacx.swcx.frozen = True
m.vacx.model.frozen = True
m.vacx.norm.frozen = True
m.vacx_12.kT.frozen = True
m.vacx_12.FracHe0.link = "72"
m.vacx_12.redshift.frozen = True
m.vacx_12.swcx.frozen = True
m.vacx_12.model.frozen = True
m.vacx_12.norm.frozen = True
m.vacx_13.kT.frozen = True
m.vacx_13.FracHe0.link = "72"
m.vacx_13.redshift.frozen = True
m.vacx_13.swcx.frozen = True
m.vacx_13.model.frozen = True
m.vacx_13.norm.frozen = True
m.vacx_14.kT.frozen = True
m.vacx_14.FracHe0.link = "72"
m.vacx_14.redshift.frozen = True
m.vacx_14.swcx.frozen = True
m.vacx_14.model.frozen = True
m.vacx_14.norm.frozen = True

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

# thaw CX
m.vacx.kT.frozen = False
m.vacx.norm.frozen = False
m.vacx_12.kT.frozen = False
m.vacx_12.norm.frozen = False
m.vacx_13.kT.frozen = False
m.vacx_13.norm.frozen = False
m.vacx_14.kT.frozen = False
m.vacx_14.norm.frozen = False

xs.Fit.perform()

# set plotting parameters
xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")

energy = 1000*np.array(xs.Plot.x())
counts = np.array(xs.Plot.y())
model = np.array(xs.Plot.model())
mpl.fill_between(energy,counts,step='mid',alpha=.4)
mpl.step(energy,counts,where="mid",lw=.5)
mpl.plot(energy,model,lw=2)
mpl.ylabel("cts/sec/bin")
mpl.ylim(ymin=0)
mpl.grid(ls=":")
mpl.xlabel("Energy (eV)")
mpl.show(block=True)




