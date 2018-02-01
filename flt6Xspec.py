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

# scheme for estimating variance (needed for empty bins)
xs.Fit.weight = "model"
#xs.Fit.weight = "churazov"

# set fit statistic
xs.Fit.statMethod = "cstat"

# load spectrum
s = xs.Spectrum("k8r61.pha")

## solar abundance (Anders and Grevesse 1989)
#xs.Xset.abund = "angr"
#z = "Anders and Grevesse"
#
#m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp + apec)")
#m.setPars(3.318,.01,25,
#          3.589,.01,5,
#          .086,1,0,0,
#          1.,
#          1.54,1.2,1.4,5.7,
#          1.96,1.2,1.4,4.9,
#          .17,1,0,0)
#          
## freeze power laws
#m.bknpower.PhoIndx1.frozen = True
#m.bknpower.BreakE.frozen = True
#m.bknpower.PhoIndx2.frozen = True
##m.bknpower.norm.frozen = True
#m.bknpower_6.PhoIndx1.frozen = True
#m.bknpower_6.BreakE.frozen = True
#m.bknpower_6.PhoIndx2.frozen = True
##m.bknpower_6.norm.frozen = True
#m.bknpower_6.norm.link = "4.9 * 15 / 5.7"
#
## freeze absorbtion
#m.avgabs.norm.frozen = True
#
## freeze thermal norms
#m.apec.norm.frozen = True
#m.apec_7.norm.frozen = True
#
## freeze temperatures so they don't get out of wack
#m.apec.kT.frozen = True
#m.apec_7.kT.frozen = True
#
## set gaussian widths equal
#m.gaussian_2.Sigma.link = "2"
#
## ignore low energy bins
#s.notice("all")
#s.ignore("**-1.5 5.-**")
#
#xs.Fit.perform()
#
## freeze cal source
#m.gaussian.LineE.frozen = True
#m.gaussian.Sigma.frozen = True
#m.gaussian.norm.frozen = True
#m.gaussian_2.LineE.frozen = True
#m.gaussian_2.Sigma.frozen = True
#m.gaussian_2.norm.frozen = True
#
## freeze AGN norm
#m.bknpower.norm.frozen = True
#
## un-freeze thermal components
#m.apec.norm.frozen = False
#m.apec.kT.frozen = False
#m.apec_7.kT.frozen = False
#m.apec_7.norm.frozen = False
#
#s.notice("all")
#s.ignore("1.-**")
#
#xs.Fit.perform()
#
## access important numbers
#kb = 8.617e-8 # keV/K
#cmPerPc = 3.086e+18 # pc to cm conversion
#LHBT = m.apec.kT.values[0]/kb
#LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#HaloT = m.apec_7.kT.values[0]/kb
#HaloEM = m.apec_7.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#label = "LHB T: %.1e K\nLHB EM: %.1e pc cm^-6" % (LHBT,LHBEM)
#label += "\nHalo T: %.1e K\nHalo EM: %.1e pc cm^-6\n" % (HaloT,HaloEM)
#redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
#label += r"$X^{2}_{red}$ = %.2f" % redChi2
#
## set plotting parameters
##s.notice("all")
#xs.Plot.xAxis = "keV"
#xs.Plot.xLog = False
#xs.Plot("data")
#
#energy = 1000*np.array(xs.Plot.x())
#counts = np.array(xs.Plot.y())
#model = np.array(xs.Plot.model())
#residual = counts - model
#ylim = 1.1*np.max(np.abs(residual))
#mpl.figure()
#mpl.subplots_adjust(wspace=0,hspace=0)
#mpl.subplot(2,1,1)
#mpl.step(energy,counts,where="mid",lw=.5)
#mpl.plot(energy,model,label=label)
#mpl.ylabel("cts/sec/bin")
#mpl.ylim(ymin=0)
##mpl.yscale('log')
#mpl.title("APEC Halo\n%s Abundances" % z)
#mpl.legend(loc='best')
#mpl.grid(ls=":")
#mpl.subplot(4,1,3)
#mpl.xlabel("Energy (eV)")
#mpl.step(energy,residual,where="mid",lw=.5)
#mpl.grid(ls=":")
#mpl.ylim(-1*ylim,ylim)
#mpl.axhline(0)
#mpl.show(block=True)
#
#################################################################################
#
## cool cloud depleted abundances (Savage and Sembach 1996)
#xs.Xset.abund = "file sscooln.abund"
#z = "Savage & Sembach"
#
#m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp + apec)")
#m.setPars(3.318,.01,25,
#          3.589,.01,5,
#          .086,1,0,0,
#          1.,
#          1.54,1.2,1.4,5.7,
#          1.96,1.2,1.4,4.9,
#          .17,1,0,0)
#          
## freeze power laws
#m.bknpower.PhoIndx1.frozen = True
#m.bknpower.BreakE.frozen = True
#m.bknpower.PhoIndx2.frozen = True
##m.bknpower.norm.frozen = True
#m.bknpower_6.PhoIndx1.frozen = True
#m.bknpower_6.BreakE.frozen = True
#m.bknpower_6.PhoIndx2.frozen = True
##m.bknpower_6.norm.frozen = True
#m.bknpower_6.norm.link = "4.9 * 15 / 5.7"
#
## freeze absorbtion
#m.avgabs.norm.frozen = True
#
## freeze thermal norms
#m.apec.norm.frozen = True
#m.apec_7.norm.frozen = True
#
## freeze temperatures so they don't get out of wack
#m.apec.kT.frozen = True
#m.apec_7.kT.frozen = True
#
## set gaussian widths equal
#m.gaussian_2.Sigma.link = "2"
#
## ignore low energy bins
#s.notice("all")
#s.ignore("**-1.5 5.-**")
#
#xs.Fit.perform()
#
## freeze cal source
#m.gaussian.LineE.frozen = True
#m.gaussian.Sigma.frozen = True
#m.gaussian.norm.frozen = True
#m.gaussian_2.LineE.frozen = True
#m.gaussian_2.Sigma.frozen = True
#m.gaussian_2.norm.frozen = True
#
## freeze AGN norm
#m.bknpower.norm.frozen = True
#
## un-freeze thermal components
#m.apec.norm.frozen = False
#m.apec.kT.frozen = False
#m.apec_7.kT.frozen = False
#m.apec_7.norm.frozen = False
#
#s.notice("all")
#s.ignore("1.-**")
#
#xs.Fit.perform()
#
## access important numbers
#kb = 8.617e-8 # keV/K
#cmPerPc = 3.086e+18 # pc to cm conversion
#LHBT = m.apec.kT.values[0]/kb
#LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#HaloT = m.apec_7.kT.values[0]/kb
#HaloEM = m.apec_7.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#label = "LHB T: %.1e K\nLHB EM: %.1e pc cm^-6" % (LHBT,LHBEM)
#label += "\nHalo T: %.1e K\nHalo EM: %.1e pc cm^-6\n" % (HaloT,HaloEM)
#redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
#label += r"$X^{2}_{red}$ = %.2f" % redChi2
#
## set plotting parameters
##s.notice("all")
#xs.Plot.xAxis = "keV"
#xs.Plot.xLog = False
#xs.Plot("data")
#
#energy = 1000*np.array(xs.Plot.x())
#counts = np.array(xs.Plot.y())
#model = np.array(xs.Plot.model())
#residual = counts - model
#ylim = 1.1*np.max(np.abs(residual))
#mpl.figure()
#mpl.subplots_adjust(wspace=0,hspace=0)
#mpl.subplot(2,1,1)
#mpl.step(energy,counts,where="mid",lw=.5)
#mpl.plot(energy,model,label=label)
#mpl.ylabel("cts/sec/bin")
#mpl.ylim(ymin=0)
##mpl.yscale('log')
#mpl.title("APEC Halo\n%s Abundances" % z)
#mpl.legend(loc='best')
#mpl.grid(ls=":")
#mpl.subplot(4,1,3)
#mpl.xlabel("Energy (eV)")
#mpl.step(energy,residual,where="mid",lw=.5)
#mpl.grid(ls=":")
#mpl.ylim(-1*ylim,ylim)
#mpl.axhline(0)
#mpl.show(block=True)
#
#################################################################################
#
#xs.Xset.abund = "angr"
#z = "Fujimoto 2007"
#
#m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp + vapec)")
#m.setPars(3.318,.01,25,
#          3.589,.01,5,
#          .086,1,0,0,
#          1.,
#          1.54,1.2,1.4,5.7,
#          1.96,1.2,1.4,4.9,
#          .17,1,1.92,2.14,1,2.77,1,1,1,1,1,1,1.42,1,0,0,)
#          
## freeze power laws
#m.bknpower.PhoIndx1.frozen = True
#m.bknpower.BreakE.frozen = True
#m.bknpower.PhoIndx2.frozen = True
##m.bknpower.norm.frozen = True
#m.bknpower_6.PhoIndx1.frozen = True
#m.bknpower_6.BreakE.frozen = True
#m.bknpower_6.PhoIndx2.frozen = True
##m.bknpower_6.norm.frozen = True
#m.bknpower_6.norm.link = "4.9 * 15 / 5.7"
#
## freeze absorbtion
#m.avgabs.norm.frozen = True
#
## freeze thermal norms
#m.apec.norm.frozen = True
#m.vapec.norm.frozen = True
#
## freeze temperatures so they don't get out of wack
#m.apec.kT.frozen = True
#m.vapec.kT.frozen = True
#
## set gaussian widths equal
#m.gaussian_2.Sigma.link = "2"
#
## ignore low energy bins
#s.notice("all")
#s.ignore("**-1.5 5.-**")
#
#xs.Fit.perform()
#
## freeze cal source
#m.gaussian.LineE.frozen = True
#m.gaussian.Sigma.frozen = True
#m.gaussian.norm.frozen = True
#m.gaussian_2.LineE.frozen = True
#m.gaussian_2.Sigma.frozen = True
#m.gaussian_2.norm.frozen = True
#
## freeze AGN norm
#m.bknpower.norm.frozen = True
#
## un-freeze thermal components
#m.apec.norm.frozen = False
#m.apec.kT.frozen = False
#m.vapec.kT.frozen = False
#m.vapec.norm.frozen = False
#
#s.notice("all")
#s.ignore("1.-**")
#
#xs.Fit.perform()
#
## access important numbers
#kb = 8.617e-8 # keV/K
#cmPerPc = 3.086e+18 # pc to cm conversion
#LHBT = m.apec.kT.values[0]/kb
#LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#HaloT = m.vapec.kT.values[0]/kb
#HaloEM = m.vapec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#label = "LHB T: %.1e K\nLHB EM: %.1e pc cm^-6" % (LHBT,LHBEM)
#label += "\nHalo T: %.1e K\nHalo EM: %.1e pc cm^-6\n" % (HaloT,HaloEM)
#redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
#label += r"$X^{2}_{red}$ = %.2f" % redChi2
#
## set plotting parameters
##s.notice("all")
#xs.Plot.xAxis = "keV"
#xs.Plot.xLog = False
#xs.Plot("data")
#
#energy = 1000*np.array(xs.Plot.x())
#counts = np.array(xs.Plot.y())
#model = np.array(xs.Plot.model())
#residual = counts - model
#ylim = 1.1*np.max(np.abs(residual))
#mpl.figure()
#mpl.subplots_adjust(wspace=0,hspace=0)
#mpl.subplot(2,1,1)
#mpl.step(energy,counts,where="mid",lw=.5)
#mpl.plot(energy,model,label=label)
#mpl.ylabel("cts/sec/bin")
#mpl.ylim(ymin=0)
##mpl.yscale('log')
#mpl.title("APEC Halo\n%s Abundances" % z)
#mpl.legend(loc='best')
#mpl.grid(ls=":")
#mpl.subplot(4,1,3)
#mpl.xlabel("Energy (eV)")
#mpl.step(energy,residual,where="mid",lw=.5)
#mpl.grid(ls=":")
#mpl.ylim(-1*ylim,ylim)
#mpl.axhline(0)
#mpl.show(block=True)
################################################################################
#
#xs.Xset.abund = "angr"
#
#m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp) + acx")
#m.setPars(3.318,.01,25,
#          3.589,.01,5,
#          .08617,1,0,0,
#          1.,
#          1.54,1.2,1.4,5.7,
#          1.96,1.2,1.4,4.9,
#          .17,.1,1,0,1,8,0)
#          
## freeze power laws
#m.bknpower.PhoIndx1.frozen = True
#m.bknpower.BreakE.frozen = True
#m.bknpower.PhoIndx2.frozen = True
##m.bknpower.norm.frozen = True
#m.bknpower_6.PhoIndx1.frozen = True
#m.bknpower_6.BreakE.frozen = True
#m.bknpower_6.PhoIndx2.frozen = True
##m.bknpower_6.norm.frozen = True
#m.bknpower_6.norm.link = "4.9 * 15 / 5.7"
#
## freeze absorbtion
#m.avgabs.norm.frozen = True
#
## freeze thermal norms
#m.apec.norm.frozen = True
#
## freeze temperatures so they don't get out of wack
#m.apec.kT.frozen = True
#
## freeze CX
#m.acx.kT.frozen = True
#m.acx.FracHe0.frozen = True
#m.acx.Abundanc.frozen = True
#m.acx.redshift.frozen = True
#m.acx.swcx.frozen = True
#m.acx.model.frozen = True
#m.acx.norm.frozen = True
#
## set gaussian widths equal
#m.gaussian_2.Sigma.link = "2"
#
## ignore low energy bins
#s.notice("all")
#s.ignore("**-1.5 5.-**")
#
#xs.Fit.perform()
#
## freeze cal source
#m.gaussian.LineE.frozen = True
#m.gaussian.Sigma.frozen = True
#m.gaussian.norm.frozen = True
#m.gaussian_2.LineE.frozen = True
#m.gaussian_2.Sigma.frozen = True
#m.gaussian_2.norm.frozen = True
#
## freeze AGN norm
#m.bknpower.norm.frozen = True
#
## un-freeze thermal components
#m.apec.norm.frozen = False
#m.apec.kT.frozen = False
#
## un-freeze CX
#m.acx.kT.frozen = False
#m.acx.FracHe0.frozen = False
#m.acx.norm.frozen = False
#
#s.notice("all")
#s.ignore("1.-**")
#
#xs.Fit.perform()
#
## access important numbers
#kb = 8.617e-8 # keV/K
#cmPerPc = 3.086e+18 # pc to cm conversion
#LHBT = m.apec.kT.values[0]/kb
#LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#label = "LHB T: %.1e* K\nLHB EM: %.1e pc cm^-6\n" % (LHBT,LHBEM)
#redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
#label += r"$X^{2}_{red}$ = %.2f" % redChi2
#
## set plotting parameters
#xs.Plot.xAxis = "keV"
#xs.Plot.xLog = False
#xs.Plot("data")
#
#energy = 1000*np.array(xs.Plot.x())
#counts = np.array(xs.Plot.y())
#model = np.array(xs.Plot.model())
#residual = counts - model
#ylim = 1.1*np.max(np.abs(residual))
#mpl.figure()
#mpl.subplots_adjust(wspace=0,hspace=0)
#mpl.subplot(2,1,1)
#mpl.step(energy,counts,where="mid",lw=.5)
#mpl.plot(energy,model,label=label)
#mpl.ylabel("cts/sec/bin")
#mpl.ylim(ymin=0)
#mpl.title("No Halo\nSWCX Only")
#mpl.legend(loc='best')
#mpl.grid(ls=":")
#mpl.subplot(4,1,3)
#mpl.xlabel("Energy (eV)")
#mpl.step(energy,residual,where="mid",lw=.5)
#mpl.grid(ls=":")
#mpl.ylim(-1*ylim,ylim)
#mpl.axhline(0)
#mpl.show(block=True)
#
###############################################################################


xs.Xset.abund = "angr"

m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp + apec) + vacx + vacx + vacx + vacx")
m.setPars(3.318,.01,25,
          3.589,.01,5,
          .08617,1,0,0,
          1.,
          1.54,1.2,1.4,5.7,
          1.96,1.2,1.4,4.9,
          .17,1,0,0,
          .08,.1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,8,0,
          .17,.1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,8,0,
          1.,.1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,8,0,
          .17,.1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,8,0)
          
# freeze power laws
m.bknpower.PhoIndx1.frozen = True
m.bknpower.BreakE.frozen = True
m.bknpower.PhoIndx2.frozen = True
#m.bknpower.norm.frozen = True
m.bknpower_6.PhoIndx1.frozen = True
m.bknpower_6.BreakE.frozen = True
m.bknpower_6.PhoIndx2.frozen = True
#m.bknpower_6.norm.frozen = True
m.bknpower_6.norm.link = "4.9 * 15 / 5.7"

# freeze absorbtion
m.avgabs.norm.frozen = True

# freeze thermal norms
m.apec.norm.frozen = True
m.apec_7.norm.frozen = True

# freeze temperatures so they don't get out of wack
m.apec.kT.frozen = True
m.apec_7.kT.frozen = True

# freeze CX
m.vacx.kT.frozen = True
m.vacx.FracHe0.frozen = True
m.vacx.redshift.frozen = True
m.vacx.swcx.frozen = True
m.vacx.model.frozen = True
m.vacx.norm.frozen = True
m.vacx_9.kT.frozen = True
m.vacx_9.FracHe0.link = "25"
m.vacx_9.redshift.frozen = True
m.vacx_9.swcx.frozen = True
m.vacx_9.model.frozen = True
m.vacx_9.norm.frozen = True
m.vacx_10.kT.frozen = True
m.vacx_10.FracHe0.link = "25"
m.vacx_10.redshift.frozen = True
m.vacx_10.swcx.frozen = True
m.vacx_10.model.frozen = True
m.vacx_10.norm.frozen = True
m.vacx_11.kT.frozen = True
m.vacx_11.FracHe0.link = "25"
m.vacx_11.redshift.frozen = True
m.vacx_11.swcx.frozen = True
m.vacx_11.model.frozen = True
m.vacx_11.norm.frozen = True

# set gaussian widths equal
m.gaussian_2.Sigma.link = "2"

# ignore low energy bins
s.notice("all")
s.ignore("**-1.5 5.-**")

xs.Fit.perform()

# freeze cal source
m.gaussian.LineE.frozen = True
m.gaussian.Sigma.frozen = True
m.gaussian.norm.frozen = True
m.gaussian_2.LineE.frozen = True
m.gaussian_2.Sigma.frozen = True
m.gaussian_2.norm.frozen = True

# freeze AGN norm
m.bknpower.norm.frozen = True

# un-freeze thermal components
m.apec.norm.frozen = False
m.apec.kT.frozen = False
m.apec_7.norm.frozen = False
m.apec_7.kT.frozen = False

s.notice("all")
s.ignore("**-.2 1.-**")

xs.Fit.perform()

# access important numbers
kb = 8.617e-8 # keV/K
cmPerPc = 3.086e+18 # pc to cm conversion
LHBT = m.apec.kT.values[0]/kb
LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
label = r"LHB T: %.1e K" % LHBT
label += "\n" + r"LHB EM: %.1e $cm^{-6}pc$" % LHBEM
HaloT = m.apec_7.kT.values[0]/kb
HaloEM = m.apec_7.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
label += "\n" + r"Halo T: %.1e K" % HaloT
label += "\n" + "Halo EM: %.1e $cm^{-6}pc$" % HaloEM
redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
label += "\n" + r"$X^{2}_{red}$ = %.2f" % redChi2

# set plotting parameters
xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")

energy = 1000*np.array(xs.Plot.x())
counts = np.array(xs.Plot.y())
model = np.array(xs.Plot.model())
residual = counts - model
ylim = 1.1*np.max(np.abs(residual))
mpl.figure()
mpl.subplots_adjust(wspace=0,hspace=0)
mpl.subplot(2,1,1)
mpl.step(energy,counts,where="mid",lw=.5)
mpl.plot(energy,model,label=label)
mpl.ylabel("cts/sec/bin")
mpl.ylim(ymin=0)
mpl.title("Halo Only")
mpl.legend(loc='best')
mpl.grid(ls=":")
mpl.subplot(4,1,3)
mpl.xlabel("Energy (eV)")
mpl.step(energy,residual,where="mid",lw=.5)
mpl.grid(ls=":")
mpl.ylim(-1*ylim,ylim)
mpl.axhline(0)
mpl.show(block=True)

# un-freeze CX
m.vacx.kT.frozen = False
m.vacx.FracHe0.frozen = False
m.vacx.norm.frozen = False
m.vacx_9.kT.frozen = False
m.vacx_9.norm.frozen = False
m.vacx_10.kT.frozen = False
m.vacx_10.norm.frozen = False
m.vacx_11.kT.frozen = False
m.vacx_11.norm.frozen = False

s.notice("all")
s.ignore("**-.2 1.5-**")

xs.Fit.perform()

# access important numbers
kb = 8.617e-8 # keV/K
cmPerPc = 3.086e+18 # pc to cm conversion
LHBT = m.apec.kT.values[0]/kb
LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
label = r"LHB T: %.1e K" % LHBT
label += "\n" + r"LHB EM: %.1e $cm^{-6}pc$" % LHBEM
HaloT = m.apec_7.kT.values[0]/kb
HaloEM = m.apec_7.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
label += "\n" + r"Halo T: %.1e K" % HaloT
label += "\n" + "Halo EM: %.1e $cm^{-6}pc$" % HaloEM
redChi2 =  xs.Fit.testStatistic/xs.Fit.dof
label += "\n" + r"$X^{2}_{red}$ = %.2f" % redChi2

# set plotting parameters
xs.Plot.xAxis = "keV"
xs.Plot.xLog = False
xs.Plot("data")

energy = 1000*np.array(xs.Plot.x())
counts = np.array(xs.Plot.y())
model = np.array(xs.Plot.model())
residual = counts - model
ylim = 1.1*np.max(np.abs(residual))
mpl.figure()
mpl.subplots_adjust(wspace=0,hspace=0)
mpl.subplot(2,1,1)
mpl.step(energy,counts,where="mid",lw=.5)
mpl.plot(energy,model,label=label)
mpl.ylabel("cts/sec/bin")
mpl.ylim(ymin=0)
mpl.title("Halo + SWCX")
mpl.legend(loc='best')
mpl.grid(ls=":")
mpl.subplot(4,1,3)
mpl.xlabel("Energy (eV)")
mpl.step(energy,residual,where="mid",lw=.5)
mpl.grid(ls=":")
mpl.ylim(-1*ylim,ylim)
mpl.axhline(0)
mpl.show(block=True)

LHBN = m.apec.norm.values[0]
HaloN = m.apec_7.norm.values[0]
AGN1N = m.bknpower.norm.values[0]
AGN2N = m.bknpower.norm.values[0]
CN = m.vacx.norm.values[0]
ON = m.vacx_9.norm.values[0]
NeN = m.vacx_10.norm.values[0]
FeN = m.vacx_11.norm.values[0]






##############################################################################
#
#xs.Xset.abund = "file sscooln.abund"
#
#m = xs.Model("gau + gau + apec + avgabs*(bknp + bknp + apec) + acx")
#m.setPars(3.318,.01,25,
#          3.589,.01,5,
#          .103192,1,0,5.25123,
#          1.,
#          1.54,1.2,1.4,5.7,
#          1.96,1.2,1.4,4.9,
#          .239310,1,0,17.7314,
#          .258162,.252242,1,0,1,8,.0385316)
#          
## freeze power laws
#m.bknpower.PhoIndx1.frozen = True
#m.bknpower.BreakE.frozen = True
#m.bknpower.PhoIndx2.frozen = True
#m.bknpower.norm.frozen = True
#m.bknpower_6.PhoIndx1.frozen = True
#m.bknpower_6.BreakE.frozen = True
#m.bknpower_6.PhoIndx2.frozen = True
#m.bknpower_6.norm.frozen = True
#
## freeze absorbtion
#m.avgabs.norm.frozen = True
#
## freeze thermal norms
#m.apec.norm.frozen = True
#m.apec_7.norm.frozen = True
#
## freeze temperatures so they don't get out of wack
#m.apec.kT.frozen = True
#m.apec_7.kT.frozen = True
#
## freeze CX
#m.acx.kT.frozen = True
#m.acx.FracHe0.frozen = True
#m.acx.Abundanc.frozen = True
#m.acx.redshift.frozen = True
#m.acx.swcx.frozen = True
#m.acx.model.frozen = True
#m.acx.norm.frozen = True
#
## set gaussian widths equal
#m.gaussian_2.Sigma.link = "2"
#
## ignore low energy bins
#s.notice("all")
#s.ignore("**-3.2")
#
#xs.Fit.perform()
#
## freeze cal source
#m.gaussian.LineE.frozen = True
#m.gaussian.Sigma.frozen = True
#m.gaussian.norm.frozen = True
#m.gaussian_2.LineE.frozen = True
#m.gaussian_2.Sigma.frozen = True
#m.gaussian_2.norm.frozen = True
#
#s.notice("all")
#s.ignore("1.-**")
#
#m.apec.norm.frozen = False
#m.apec_7.norm.frozen = False
##m.apec.kT.frozen = False
#m.apec_7.kT.frozen = False
##m.acx.kT.frozen = False
##m.acx.FracHe0.frozen = False
#m.acx.norm.frozen = False
#
## step through fixed ratios of halo and SWCX
#xs.Fit.steppar("23 0. 17.7314 10 30 0. .03853 10")
#xs.Plot.device = "/xs"
#xs.Plot("contour")

#xs.Plot("data")
#
#s.notice("all")
#s.ignore("1.2-**")
#
#xs.Fit.perform()
#
## Interpretation: SWCX-only not favored, 
#
## access important numbers
#kb = 8.617e-8 # keV/K
#cmPerPc = 3.086e+18 # pc to cm conversion
#LHBT = m.apec.kT.values[0]/kb
#LHBEM = m.apec.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#HaloT = m.apec_7.kT.values[0]/kb
#HaloEM = m.apec_7.norm.values[0]*4*np.pi/(1e-14*cmPerPc)
#label = "LHB T: %.1e* K\nLHB EM: %.1e pc cm^-6" % (LHBT,LHBEM)
#label += "\nHalo T: %.1e K\nHalo EM: %.1e pc cm^-6" % (HaloT,HaloEM)
#
## set plotting parameters
#xs.Plot.xAxis = "keV"
#xs.Plot.xLog = False
#xs.Plot("data")
#
#energy = 1000*np.array(xs.Plot.x())
#counts = np.array(xs.Plot.y())
#model = np.array(xs.Plot.model())
#residual = counts - model
#ylim = 1.1*np.max(np.abs(residual))
#mpl.figure()
#mpl.subplots_adjust(wspace=0,hspace=0)
#mpl.subplot(2,1,1)
#mpl.step(energy,counts,where="mid",lw=.5)
#mpl.plot(energy,model,label=label)
#mpl.ylabel("cts/sec/bin")
#mpl.ylim(ymin=0)
#mpl.title("Halo + SWCX")
#mpl.legend(loc='best')
#mpl.grid(ls=":")
#mpl.subplot(4,1,3)
#mpl.xlabel("Energy (eV)")
#mpl.step(energy,residual,where="mid",lw=.5)
#mpl.grid(ls=":")
#mpl.ylim(-1*ylim,ylim)
#mpl.axhline(0)
#mpl.show(block=True)
#
#
#
