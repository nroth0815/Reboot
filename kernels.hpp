#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

MyFloat alpha(int q1, int q2, int q3, int p1, int p2, int p3){
	
	int modp=p1*p1+p2*p2+p3*p3;

return modp==0 ? 0. : (MyFloat)(q1*p1+q2*p2+q3*p3)/modp; //returns 0 if modp=0, qp/modp otherwise
}

MyFloat beta(int q1, int q2, int q3, int p1, int p2, int p3){

	int modq, modp, qp;
	MyFloat value;

	modq = q1*q1 + q2*q2 + q3*q3;
	modp = p1*p1 + p2*p2 + p3*p3;
	MyFloat invq;
	MyFloat invp;		

	if(modq == 0 || modp == 0){
			value=0.0;
	}
	else{
		qp=q1*p1+q2*p2+q3*p3;
		invq = 1.0/modq;
		invp = 1.0/modp;
		value= qp * (invq + invp + 2.0*qp*invq*invp)/2.0;
	}
	
	return value;
}

