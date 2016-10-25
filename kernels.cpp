#include "header.hpp"

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

MyFloat kernel(int q1, int q2, int q3, int p1, int p2, int p3){ //this takes k1, k2, not k=k1+k2
	//this is a non-optimized (yet correct) version for testing purposes only

		MyFloat eps, value=0.,modq,modp,qp;

#ifdef DOUBLEPRECISION
eps=1e-15;
#else
eps=1e-07;
#endif

		modq=(MyFloat)(q1*q1+q2*q2+q3*q3);
		modp=(MyFloat)(p1*p1+p2*p2+p3*p3);

		qp=(MyFloat)(q1*p1+q2*p2+q3*p3);

		value=5./7+2./7*qp*qp/(modq+eps)/(modp+eps)+1./2*qp*(1./(modq+eps)+1./(modp+eps)) ;

	return value; 
}
