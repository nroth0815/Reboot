#ifndef HEADER_H
#define HEADER_H

#include <cmath>
#include <iostream> //for cerr
#include <iomanip>
#include <fstream>

#ifdef MAC
#  include <mach/error.h>
#  include "/Users/nroth/Projects/code/HDF_IO.hh"
#  ifdef DOUBLEPRECISION
#    include <drfftw.h>
#  else
#    include <srfftw.h>
#  endif
#else
#include <error.h>
#include <rfftw.h> //probably double precision; include switch here too (need to link to localcode)
//#include "HDF_IO.hh" //curently broken somehow
#endif

#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

MyFloat alpha(int q1, int q2, int q3, int p1, int p2, int p3);
MyFloat beta(int q1, int q2, int q3, int p1, int p2, int p3);
MyFloat kernel(int q1, int q2, int q3, int p1, int p2, int p3);

int Ps(int res, MyFloat Boxlength, int nBins, fftw_complex *ft, const char *out);
fftw_complex *Smooth(const char *Outputfile, const char *Outputfile2, int res, 					 fftw_complex *datain, MyFloat R, MyFloat Boxlength);

#endif