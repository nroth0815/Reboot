#include "stats.hpp"

int Ps(int res, MyFloat Boxlength, int nBins, fftw_complex *ft, const char *out){

	MyFloat *inBin = new MyFloat[nBins];
	MyFloat *Gx = new MyFloat[nBins];
	MyFloat *kbin = new MyFloat[nBins];
	MyFloat kmax = M_PI/Boxlength*(MyFloat)res, kmin = 2.0f*M_PI/(MyFloat)Boxlength, dklog = log10(kmax/kmin)/nBins, kw = 2.0f*M_PI/(MyFloat)Boxlength;
	int ix,iy,iz;
	MyFloat kfft;
	size_t idx2,idx;

	for( ix=0; ix<nBins; ix++ ){
		inBin[ix] = 0;
		Gx[ix] = 0.0f;
		kbin[ix] = 0.0f;
	}

	int iix, iiy, iiz, r2=res/2;

	for(ix=0; ix<res;ix++){
		for(iy=0;iy<res;iy++){
			for(iz=0;iz<res;iz++){
				idx = (ix*res+iy)*(res)+iz;

				// determine mode modulus
				MyFloat vabs = ft[idx].re*ft[idx].re+ft[idx].im*ft[idx].im; //factor 2??

				if( ix>r2 ) iix = ix - res; else iix = ix;
				if( iy>r2 ) iiy = iy - res; else iiy = iy;
				if( iz>r2 ) iiz = iz - res; else iiz = iz;
					//iiz = iz;
				
				kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);
				MyFloat k = kfft*kw;
				
					// correct for aliasing, formula from Jing (2005), ApJ 620, 559
					// assume isotropic aliasing (approx. true for k<kmax=knyquist)
					// this formula is for CIC interpolation scheme, which we use <- only needed for Powerspectrum??
				MyFloat JingCorr = (1.0f-2.0f/3.0f*sin(M_PI*k/kmax/2.0f)*sin(M_PI*k/kmax/2.0f));	
				vabs /= JingCorr;

					//.. logarithmic spacing in k
				idx2 = (size_t)(( 1.0f/dklog * log10(k/kmin) ));
				
				if(k>=kmin&&k<kmax){
					Gx[idx2] += vabs;
					kbin[idx2] += k;
					inBin[idx2]++;
				}	
			}
		}
	}
   
	MyFloat psnorm2=pow((2.0*M_PI)/Boxlength,3.0); //whatever, don't ask, this works
	std::ofstream ofs(out);
	for(ix=0; ix<nBins; ix++){
		if(inBin[ix]>0){
			ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) )
			<< std::setw(16) <<  kbin[ix]/inBin[ix]
			<< std::setw(16) << (MyFloat)(Gx[ix]/inBin[ix])/psnorm2
			<< std::setw(16) << inBin[ix]
			<< std::endl;

		}
		else{continue;}
	}

   	ofs.close();
   	delete[] inBin;
   	delete[] Gx;
   	delete[] kbin;

   	return 0;
}

fftw_complex* Smooth(const char *Outputfile, const char *Outputfile2, int res, fftw_complex *datain, MyFloat R, MyFloat Boxlength){

	fftw_complex *data = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	fftw_complex *deltaw=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));

	
	int i, k1, k2, k3, kk1,kk2,kk3, r2=res/2;

  	MyFloat k,w,kw = 2.0f*M_PI/(MyFloat)Boxlength; //kmin
  	std::ofstream out2(Outputfile2);

	//calculate mod(k) and multiply with smoothing function

    for(k1=0; k1<res;k1++)
	  for(k2=0;k2<res;k2++)
		for(k3=0;k3<res;k3++){
			i = (k1*res+k2)*(res)+k3;

			if( k1>r2 ) kk1 = k1-res; else kk1 = k1;
			if( k2>r2 ) kk2 = k2-res; else kk2 = k2;
			if( k3>r2 ) kk3 = k3-res; else kk3 = k3;		

			k=(MyFloat)(kk1*kk1+kk2*kk2+kk3*kk3)*kw*kw;
			w=(MyFloat)exp((MyFloat)(-(MyFloat)(k*R*R)/2.));

			deltaw[i].re=w*datain[i].re;	
			deltaw[i].im=w*datain[i].im; 
			out2 << deltaw[i].re<<" "<< deltaw[i].im << std::endl;
	}

	out2.close();

	//fftw smoothed data back to real-space
	fftwnd_plan p;
	p = fftw3d_create_plan( res,res,res, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftwnd_one(p, deltaw, data);
	fftwnd_destroy_plan(p);


	std::ofstream out(Outputfile);
	for(i=0;i<res*res*res;i++){ out << data[i].re << std::endl;}
	out.close();

	fftw_free(data);

	return deltaw;
}