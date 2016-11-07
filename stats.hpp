#ifndef STATS_H
#define STATS_H

#include <cmath>
#include <iomanip> //for setw
#include <fstream>

#ifdef MAC
#  include "/Users/nroth/Projects/code/HDF_IO.hh"
#  ifdef DOUBLEPRECISION
#    include <drfftw.h>
#  else
#    include <srfftw.h>
#  endif
#else
#include <rfftw.h> //probably double precision; include switch here too (need to link to localcode)
//#include "HDF_IO.hh" //curently broken somehow
#endif

#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

int Ps(int res, MyFloat Boxlength, int nBins, fftw_complex *ft, const char *out);
fftw_complex *Smooth(const char *Outputfile, const char *Outputfile2, int res,
 					 fftw_complex *datain, MyFloat R, MyFloat Boxlength);

#endif