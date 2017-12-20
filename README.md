# flights
By Dallas Wulf

This is the directory where I have tried to contain all of 
the data and code needed to analyze XQC flight data.

flightLibrary.py 
This is a modified version of general filterLibrary.py, and 
it contains all of the run-specific definitions and numbers 
for self-contained processing of flight data. This file only
contains classes and defintions that are called by other 
programs. It doesn't do anything on its own.

analyzeFlight6.py
Complete analysis of run k8r61 a.k.a. Flight 6 a.k.a. 
Flight 36.294 a.k.a. the 2013 flight.
All data dependecies go in the /data directory:
* 6 cardfiles (k8r61*.le.dat)
* cold-plate temperature (tcp_Iflight)
* optimum filter .fits file (k8O61_fil.fits)
* optimum filter .fits file for frame pixel (k8F61_fil.fits)
* 6 IR filter transmission measurements:
- 11895-2.dat
- ODB103_092013.dat
- 2K_30048-3.dat (this is not a very good measurement)
- 2K105_082013.dat 
- 130K116_082013.dat
- RT114_092013.dat
* /length directory (column depths for fitting IR filters):
- Al_smith.len
- O.len
- Parylene.len
- PI.len
- Si.len

This script packs the 6 card files into k8r61.bsn.hdf5, applies
timing from the recored bias-offs, optimally filters and fits 
all the pulses in the data stream, uses the cold plate temperature 
to apply a non-linear energy gain correction to the pulse heights,
fits the filter stack transmission, selects good events and 
calculates live time, and finally produces a spectrum and response 
file that can be read by xspec. In the future, I hope to include the 
pyxspec fitting in this script.  All processing information is stored 
in the bsn.hdf5 file, which can be easily accessed through the 
singlePixel class in flightLibrary.py

genrsp.tmplt, genpha.tmplt
These are templates for the heasarc utilties need to make PHA and 
RMF/RSP files. They also include all hte documentation I could find 
online, since it wasn't obvious to me how to use the utilties.

analyzeBelljar.py
This script is very similar to the first part of analyzeFlight6.py, 
except for l5r18 (belljar source run). However, rather than
making an xspec spectrum, this is just used to determine the 
photoelectic escape fraction.  Also, instead of using bias-offs,
this script finds sweeps in the data stream to get the timing.
/data directory dependecies:
* 6 card files (l5r18*.le.dat)
* swp file for timing (l5r18_swp.fits)

analyzeFlight5.py
Complete analysis of run j5r62 a.k.a. Flight 5 a.k.a.
Flight 36.264 a.k.a. the 2011 flight.
/data dependecies:
* 6 card files (j5r62*.le.dat)
* swp file for timing (j5r62_swp.fits)



2008 Flight = 36.223 = h6r19
