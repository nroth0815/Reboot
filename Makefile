#CC    	=  mpiCC
CC = /usr/bin/g++
#CC      =  /opt/local/bin/g++ #this is the more up-to-date and system-standard C11 compiler but it doesn't work here
#mpiCC
CFLAGS	= -Wall -O3
#-O2
#CPATH = /users/giannant/local/include
CPATH = /usr/local/include
CPATH2= /users/nroth/localcode/include
CPATH3= /users/nroth/localcode/fftw/include
LFLAGS  =
#LPATH   = /usr/lib64
LPATH   = /usr/local/lib
LPATH2   = /users/nroth/localcode/lib
LPATH3 = /users/nroth/localcode/fftw/lib
#LPATH   = /users/giannant/local/lib

# to show which libraries are linked to a program: ldd ./program; otool -L ./program on MAC

#-DMAC flag to include "error.h" on Mac OS

#all: ps22 ps13 ps11real ps13real delta2p_trunc
#all: delta2p_trunc delta2p_trunc_s0 dtrunc
#all: dtrunc_opt v2 smo dtrunc_opt3 delta2p

all: delta2pclvhm
#all: delta2p delta2pclv hm
#all: MACdelta2p MACdelta2pclv
hva: hash_v_array.cpp hash_temp.hpp Makefile
	$(CC) $(CFLAGS) -std=c++11 hash_temp.hpp -o hva hash_v_array.cpp -lm -lgsl -lgslcblas
hm: hash_test.cpp hash_temp.hpp Makefile
	$(CC) $(CFLAGS) -std=c++11 hash_temp.hpp -o hash_test hash_test.cpp -lm 
delta2p: delta2part.cpp Makefile
	$(CC) $(CFLAGS) -o delta2part -L$(LPATH2) -I$(CPATH2) kernels.hpp delta2part.cpp -lrfftw -lfftw -lm -lhdf5 -D H5_USE_16_API
delta2pclvhm: delta2part_cleverloop.cpp kernels.hpp hash_temp.hpp  Makefile
	$(CC) $(CFLAGS) -std=c++11 -o delta2part_clvhm -L$(LPATH2) -I$(CPATH2) kernels.hpp hash_temp.hpp delta2part_cleverloop.cpp -lrfftw -lfftw -lm -DDOUBLEPRECISION #-lhdf5 -D H5_USE_16_API
delta2pclv: delta2part_cleverloop.cpp kernels.hpp  Makefile
	$(CC) $(CFLAGS) -o delta2part_clv -L$(LPATH2) -I$(CPATH2) kernels.hpp delta2part_cleverloop.cpp -lrfftw -lfftw -lm -DDOUBLEPRECISION #-lhdf5 -D H5_USE_16_API
MACdelta2p: delta2part.cpp Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o delta2part -L$(LPATH) -I$(CPATH) delta2part.cpp -lsrfftw -lsfftw -lm -lhdf5 -D H5_USE_16_API -DMAC
MACdelta2pclv: delta2part_cleverloop.cpp Makefile
	$(CC) $(CFLAGS) -o delta2part_clv -L$(LPATH) -I$(CPATH) delta2part_cleverloop.cpp -lsrfftw -lsfftw -lm -lhdf5 -D H5_USE_16_API -DMAC
dptest: dp_test.cpp Makefile
	$(CC) $(CFLAGS) --std=c++0x -o dptest  -L$(LPATH) -I$(CPATH) dp_test.cpp -lm	
clean: 
	rm -rf psxx_fr psxx_half psxy_fr combsm psxx_real psxy_real
