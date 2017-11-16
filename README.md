# flights

flightLibrary.py 
Modified version of general filterLibrary.py. Contains all of 
the run-specific definitions and numbers for self-contained 
processing of flight data. Inputs are the cardfiles in 
long-endian, int16 format (i.e. the output of cardfile.f)

For k8r61, also provide fil.fits files for making optimum 
filters, and tcp_Iflight.dat for the cold-plate temperature.
(j5r62 and h6r19 make filters from the data)

For j5r62 and h6r19, also provide swp.fits files for timing.
(k8r61 uses bias offs)

All data files (cardfiles, .fits, .dat) shoould be put in data/ 
directory.

2013 Flight = 36.294 = k8r61
2011 Flight = 36.264 = j5r62
2008 Flight = 36.223 = h6r19
