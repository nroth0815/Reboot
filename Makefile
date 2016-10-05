CC    	=  /usr/bin/g++
#CC      =  /opt/local/bin/g++ #this is the more up-to-date and system-standard C11 compiler but it doesn't work here
#mpiCC
CFLAGS	= -Wall 
#-O2
#CPATH = /users/giannant/local/include
CPATH = /usr/local/include
CPATH2= /users/nroth/localcode/include
LFLAGS  =
#LPATH   = /usr/lib64
LPATH   = /usr/local/lib
LPATH2   = /users/nroth/localcode/lib
#LPATH   = /users/giannant/local/lib


#-DMAC flag to include "error.h" on Mac OS

#all: ps22 ps13 ps11real ps13real delta2p_trunc
#all: delta2p_trunc delta2p_trunc_s0 dtrunc
#all: dtrunc_opt v2 smo dtrunc_opt3 delta2p
all: dptest delta2p
delta2p: delta2part.cpp Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o delta2part -L$(LPATH) -I$(CPATH) delta2part.cpp -ldrfftw -ldfftw -lm -lhdf5 -D H5_USE_16_API -DMAC
dptest: dp_test.cpp Makefile
	$(CC) $(CFLAGS) --std=c++0x -o dptest  -L$(LPATH) -I$(CPATH) dp_test.cpp -lm 	
clean: 
	rm -rf psxx_fr psxx_half psxy_fr combsm psxx_real psxy_real
