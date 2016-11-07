#ifndef KERNELS_H
#define KERNELS_H

#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

MyFloat alpha(int q1, int q2, int q3, int p1, int p2, int p3);
MyFloat beta(int q1, int q2, int q3, int p1, int p2, int p3);
MyFloat kernel(int q1, int q2, int q3, int p1, int p2, int p3);

#endif