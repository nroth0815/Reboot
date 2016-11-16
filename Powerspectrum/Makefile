CC    	= g++
MCC= mpiCC
CFLAGS  = -Wall -O3
#CPATH = /usr/include
CPATH= /users/nroth/localcode/include
LFLAGS  =
LPATH   = /usr/lib64
LPATH2   = /users/nroth/localcode/lib
#LPATH   = /users/giannant/local/lib

all: Powerspectrum

Powerspectrum: ICBB_Powerspectrum.cc HDF_IO.hh Makefile 
	$(MCC) $(CFLAGS) $(LFLAGS) ICBB_Powerspectrum.cc -o Powerspectrum -L$(LPATH) -L$(LPATH2) -I$(CPATH) -I$(CPATH2)  -Wl,-rpath,$(LPATH2) -lsrfftw -lsfftw -Wl,-rpath,$(LPATH2) -lhdf5 -D H5_USE_16_API -lz
HaloPowerspectrum: HaloPowerspectrum.cc HDF_IO.hh Makefile
	$(MCC) $(CFLAGS) $(LFLAGS) -o HaloPowerspectrum -I$(CPATH) -L$(LPATH) HaloPowerspectrum.cc -lrfftw -lfftw -lhdf5 -lz
HaloMatterPowerspectrum: HaloMatterPowerspectrum.cc HDF_IO.hh Makefile
	$(MCC) $(CFLAGS) $(LFLAGS) -o HaloMatterPowerspectrum -I$(CPATH) -L$(LPATH) HaloMatterPowerspectrum.cc -lrfftw -lfftw -lhdf5 -lz

smooth: smoothdf.cc HDF_IO.hh Makefile 
	$(MCC) $(CFLAGS) $(LFLAGS) -o smoothdf -L$(LPATH) -I$(CPATH) smoothdf.cc -D H5_USE_16_API -lrfftw -lfftw -lhdf5 -lm

clean: 
	rm -rf *Powerspectrum
