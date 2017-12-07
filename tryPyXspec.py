# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 13:46:34 2017

@author: calorim
"""

import xspec as xs
import numpy as np
import matplotlib.pyplot as mpl

# load spectrum
s = xs.Spectrum("flt3.pha")

## define model
## cool cloud depleted abundances (Savage and Sembach 1996)
#xs.Xset.abund = "file sscooln.abund"
## unabsorbed thermal plus absorbed thermal and powerlaw
#m = xs.Model("apec + wabs*(apec + pow)")
#cmPerPc = 3.086e+18 # pc to cm conversion
#m.setPars(.107,1.,0.,.053*1e-14*cmPerPc/(4*np.pi),
#          .018,
#          .317,1.,0.,.0053*1e-14*cmPerPc/(4*np.pi),
#          1.52,12.3)
#m.show()

# solar abundance (Anders and Grevesse 1989)
xs.Xset.abund = "angr"
# unabsorbed thermal plus absorbed thermal and powerlaw
m = xs.Model("apec + wabs*(apec + pow)")
cmPerPc = 3.086e+18 # pc to cm conversion
m.setPars(.099,1.,0.,.0088*1e-14*cmPerPc/(4*np.pi),
          .018,
          .225,1.,0.,.0037*1e-14*cmPerPc/(4*np.pi),
          1.52,12.3)
m.show()

# perform fit
xs.Fit.statMethod = "cstat"
xs.Fit.perform()

# Plot with pyplot
channel = s.noticed
counts = s.values
model = m.folded(1)
mpl.figure()
mpl.step(channel, counts, where='mid', lw=.5)
mpl.plot(channel, model)
mpl.xlabel('channel')
mpl.ylabel('counts/sec/chan')
mpl.grid(ls=':')
mpl.show()