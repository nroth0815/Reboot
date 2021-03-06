#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifdef MAC 
#include <mach/error.h>
#else
#include <error.h>
#endif
#include <time.h>
using namespace std;

#ifdef MAC
#include "/Users/nroth/Projects/code/HDF_IO.hh"
#include <srfftw.h>
#else
#include <rfftw.h>
//#include "HDF_IO.hh"
#endif

int Ps(int res, int Boxlength, int nBins, fftw_complex *ft, const char *out){

	float *inBin = new float[nBins];
	float *Gx = new float[nBins];
	float *kbin = new float[nBins];
	float kmax = M_PI/Boxlength*(float)res, kmin = 2.0f*M_PI/(float)Boxlength, dklog = log10(kmax/kmin)/nBins, kw = 2.0f*M_PI/(float)Boxlength;
	int ix,iy,iz;
	float kfft;
	int idx2,idx;

	for( ix=0; ix<nBins; ix++ ){
		inBin[ix] = 0;
		Gx[ix] = 0.0f;
		kbin[ix] = 0.0f;
	}

	for(ix=0; ix<res;ix++){
		for(iy=0;iy<res;iy++){
			for(iz=0;iz<res;iz++){
				idx = (ix*res+iy)*(res)+iz;

				// determine mode modulus
				float vabs = ft[idx].re*ft[idx].re+ft[idx].im*ft[idx].im; //factor 2??

				int iix, iiy, iiz;

				if( ix>res/2 ) iix = ix - res; else iix = ix;
				if( iy>res/2 ) iiy = iy - res; else iiy = iy;
				if( iz>res/2 ) iiz = iz - res; else iiz = iz;
					//iiz = iz;
				
				kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);
				float k = kfft*kw;
				
					// correct for aliasing, formula from Jing (2005), ApJ 620, 559
					// assume isotropic aliasing (approx. true for k<kmax=knyquist)
					// this formula is for CIC interpolation scheme, which we use <- only needed for Powerspectrum??
				float JingCorr = (1.0f-2.0f/3.0f*sin(M_PI*k/kmax/2.0f)*sin(M_PI*k/kmax/2.0f));	
				vabs /= JingCorr;

					//.. logarithmic spacing in k
				idx2 = (int)(( 1.0f/dklog * log10(k/kmin) ));
				
				if(k>=kmin&&k<kmax){
					Gx[idx2] += vabs;
					kbin[idx2] += k;
					inBin[idx2]++;
				}	
			}
		}
	}
   
	float psnorm2=pow((2.0*M_PI)/Boxlength,3.0); //whatever, don't ask, this works
	ofstream ofs(out);
	for(ix=0; ix<nBins; ix++){
		if(inBin[ix]>0){
			ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) )
			<< std::setw(16) <<  kbin[ix]/inBin[ix]
			<< std::setw(16) << (float)(Gx[ix]/inBin[ix])/psnorm2
			<< std::setw(16) << inBin[ix]
			<< endl;

		}
		else{continue;}
	}

   	ofs.close();
   	delete[] inBin;
   	delete[] Gx;
   	delete[] kbin;

   	return 0;
}

fftw_complex *Smooth(const char *Outputfile, const char *Outputfile2, int res, fftw_complex *datain, float R, float Boxlength){

	fftw_complex *data = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	fftw_complex *deltaw=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));

	int k1,k2,k3,i,kk1,kk2,kk3;

  float k,w,kw = 2.0f*M_PI/(float)Boxlength; //kmin
  ofstream out2(Outputfile2);

	//calculate mod(k) and multiply with smoothing function

  for(k1=0; k1<res;k1++)
	for(k2=0;k2<res;k2++)
		for(k3=0;k3<res;k3++){
			i = (k1*res+k2)*(res)+k3;

			if( k1>res/2 ) kk1 = k1-res; else kk1 = k1;
			if( k2>res/2 ) kk2 = k2-res; else kk2 = k2;
			if( k3>res/2 ) kk3 = k3-res; else kk3 = k3;		

			k=(float)(kk1*kk1+kk2*kk2+kk3*kk3)*kw*kw;
			w=(float)exp((float)(-(float)(k*R*R)/2.));

			deltaw[i].re=w*datain[i].re;	
			deltaw[i].im=w*datain[i].im; out2<<deltaw[i].re<<" "<<deltaw[i].im<<endl;
		}

		out2.close();

	 //fftw smoothed data back to real-space
		fftwnd_plan p;
		p = fftw3d_create_plan( res,res,res, FFTW_BACKWARD, FFTW_ESTIMATE );
		fftwnd_one(p, deltaw, data);
		fftwnd_destroy_plan(p);


		ofstream out(Outputfile);
		for(i=0;i<res*res*res;i++){ out<<data[i].re<<endl;}
			out.close();

		fftw_free(data);


		return deltaw;
	}			 

float Fnew(int q1, int q2, int q3, int p1, int p2, int p3){ //this takes k1, k2, not k=k1+k2

		float eps=0.0000001,value=0.,modq,modp,qp;//, moddiff;

		modq=(float)(q1*q1+q2*q2+q3*q3);
		modp=(float)(p1*p1+p2*p2+p3*p3);
		//moddiff=(float)(diff1*diff1+diff2*diff2+diff3*diff3);

		qp=(float)(q1*p1+q2*p2+q3*p3);


		value=2./7*qp*qp/(modq+eps)/(modp+eps)+1./2*qp*(1./(modq+eps)+1./(modp+eps)) ;
		//1./2.*modq/(modp+eps)*(float)(p1*diff1+p2*diff2+p3*diff3)/(moddiff+eps);

return value; //gives correct values compared with mathematica
}

float kernel(int q1, int q2, int q3, int p1, int p2, int p3){ //this takes k1, k2, not k=k1+k2

		float eps=0.0000001,value=0.,modq,modp,qp;//, moddiff;

		modq=(float)(q1*q1+q2*q2+q3*q3);
		modp=(float)(p1*p1+p2*p2+p3*p3);
		//moddiff=(float)(diff1*diff1+diff2*diff2+diff3*diff3);

		qp=(float)(q1*p1+q2*p2+q3*p3);


		value=5./7+2./7*qp*qp/(modq+eps)/(modp+eps)+1./2*qp*(1./(modq+eps)+1./(modp+eps)) ;
		//1./2.*modq/(modp+eps)*(float)(p1*diff1+p2*diff2+p3*diff3)/(moddiff+eps);

return value; //gives correct values compared with mathematica
}

float F(int q1, int q2, int q3, int p1, int p2, int p3){

		float eps=0.0000001,value=0.,modq,modp,moddiff;

		int diff1=(q1-p1);
		int diff2=(q2-p2);
		int diff3=(q3-p3);

		modq=(float)(q1*q1+q2*q2+q3*q3);
		modp=(float)(p1*p1+p2*p2+p3*p3);
		moddiff=(float)(diff1*diff1+diff2*diff2+diff3*diff3);

		//qp=(float)(q1*p1+q2*p2+q3*p3);


		value=1./2.*modq/(modp+eps)*(float)(p1*diff1+p2*diff2+p3*diff3)/(moddiff+eps);

return value; //gives correct values compared with mathematica
}

float H(int q1, int q2, int q3, int p1, int p2, int p3){
	
	float eps=0.0000001,value=0.,modp,qp;
	
	modp=(float)(p1*p1+p2*p2+p3*p3);

	qp=(float)(q1*p1+q2*p2+q3*p3);

	value=qp/(float)(modp+eps); 


return value; //gives correct values compared with mathematica
}

int main(int argc, char *argv[]){
	
	if(argc!=7){cerr<< "Usage: ./delta2part inputfile res <smoothing scale> <'IC' or 'z0'> Boxsize <part no (0: all, or 1-8)>" <<endl; return -1;}

	string arg0=argv[0];
	string argv2=argv[2];
	string output,output2,outps;

	int res=atoi(argv[2]);
	int part=atoi(argv[6]); if(part>8 || part<0){cerr<<"part no. not between 0 and 8!"<<endl; return -1;}
	int idstop=res*res*(res/2+1),idstart=(part-1)*(idstop)/8,idstopn=idstart+(idstop)/8; if(part==0){idstart=0;idstopn=idstop;}

	int ik,jk,lk,iq1,jq1,lq1,iq2,jq2,lq2,idk,idq1,idq2,i,inq2,jnq2,lnq2,j,l,ii,jj,ll,index=0,iiq2,jjq2,llq2,id,idknew;

	int *qarr=(int*)calloc(res*res*res*3,sizeof(int));
	int *karr=(int*)calloc(idstop*4,sizeof(int));

	float R=atof(argv[3]);	
	float Boxsize=atof(argv[5]);

	float *arr=(float*)calloc(res*res*res,sizeof(float));
	float d1re, d2re, d1im, d2im,f=0.,h1=0.,h2=0.,A=0.,B=0.,C=0.;

	fftw_complex *v2=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));

	// //read in delta1(x) from gridfile and ft to delta(k)
	// const char *FILE1=argv[1];
	// int RANK=3;

	// hid_t fid1; 
	// hsize_t dims[] = {res,res,res};
	// hid_t datasetid;      // file and memory dataspace identifiers 
	// hsize_t  offset[RANK];
	
	// dims[0] = res;
	// dims[1] = res;
	// dims[2] = res;

	// offset[0] = 0;
	// offset[1] = 0;  
	// offset[2] = 0;

	// fid1 = H5Fopen(FILE1, H5F_ACC_RDONLY, H5P_DEFAULT);
	// datasetid=H5Dopen(fid1, "Density");
	// hid_t memspacebuf = H5Screate_simple(RANK, dims, NULL );
	// hid_t filespacebuf=H5Dget_space(datasetid);
	// H5Dread(datasetid,H5T_NATIVE_float, memspacebuf,filespacebuf,H5P_DEFAULT, arr);
	
	// //close identifiers
	// H5Sclose(filespacebuf);
	// H5Sclose(memspacebuf);
	// H5Dclose(datasetid);
	// H5Fclose(fid1);
	

	//normalise for fftw
	fftw_complex *arr2 = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	for(i=0;i<res*res*res;i++){arr2[i].re=arr[i]/pow(res,3.0);}//cout<<arr[i]*arr[i]<<endl;} 


	//calculate fftw
	fftw_complex *ft = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	fftwnd_plan p = fftw3d_create_plan( res, res, res, FFTW_FORWARD, FFTW_ESTIMATE );
	fftwnd_one(p, arr2, ft);
	fftwnd_destroy_plan(p);

	free(arr);
	fftw_free(arr2);

	//calulate k values
	for(i=-(res/2-1);i<res/2+1;i++){
		if(i<0){ii=res+i;}else{ii=i;}
		for(j=-(res/2-1);j<res/2+1;j++){
			if(j<0){jj=res+j;}else{jj=j;}		
			for(l=0;l<res/2+1;l++){
				
				idk=(ii*(res)+jj)*(res)+l;		
				
				karr[4*index]=i; karr[4*index+1]=j; karr[4*index+2]=l; karr[4*index+3]=idk;
				index+=1;
			}
		}
	}

	//calculate q values
	for(i=-(res/2-1);i<res/2+1;i++){
		if(i<0){ii=res+i;}else{ii=i;}
		for(j=-(res/2-1);j<res/2+1;j++){
			if(j<0){jj=res+j;}else{jj=j;}		
			for(l=-(res/2-1);l<res/2+1;l++){
				if(l<0){ll=res+l;}else{ll=l;}
				
				idk=(ii*(res)+jj)*(res)+ll;
				
				qarr[3*idk]=i; qarr[3*idk+1]=j; qarr[3*idk+2]=l; 
			}
		}
	}

	
	//float t1,t0;

	cerr<<"beginning loop"<<endl;

	clock_t t0,t1;

	//idstopn=res*res*res;//idstart+10;
	//t0= time(NULL);
	 // t0=MPI_Wtime();
	t0 = clock();

	for(id=idstart;id<idstopn;id++){
		ik=karr[id*4];
		jk=karr[id*4+1];
		lk=karr[4*id+2];
		idk=karr[4*id+3];

		
		//t0 = second();
		
		for(idq1=0;idq1<res*res*res;idq1++){
			
			iq1=qarr[idq1*3];
			jq1=qarr[idq1*3+1];
			lq1=qarr[idq1*3+2];

			//iq2=ik-iq1; if(iq2<-(res/2) ){inq2=res+iq2;} else if(iq2>res/2){inq2=iq2-res;} else{inq2=iq2;}
			//jq2=jk-jq1; if(jq2<-(res/2) ){jnq2=res+jq2;} else if(jq2>res/2){jnq2=jq2-res;} else{jnq2=jq2;}
			//lq2=lk-lq1; if(lq2<-(res/2) ){lnq2=res+lq2;} else if(lq2>res/2){lnq2=lq2-res;} else{lnq2=lq2;}

			iq2=ik-iq1; if(iq2<-(res/2) ){continue;} else if(iq2>res/2){continue;} else{inq2=iq2;}
			jq2=jk-jq1; if(jq2<-(res/2) ){continue;} else if(jq2>res/2){continue;} else{jnq2=jq2;}
			lq2=lk-lq1; if(lq2<-(res/2) ){continue;} else if(lq2>res/2){continue;} else{lnq2=lq2;}


			if(inq2<0){iiq2=res+inq2;}else{iiq2=inq2;}
			if(jnq2<0){jjq2=res+jnq2;}else{jjq2=jnq2;}
			if(lnq2<0){llq2=res+lnq2;}else{llq2=lnq2;}		

			idq2=(iiq2*res+jjq2)*res+llq2;

				//kernels
			//f=F(ik,jk,lk,iq1,jq1,lq1);
			// f=Fnew(iq1,jq1,lq1, inq2, jnq2, lnq2);
			// h1=H(ik,jk,lk,iq1,jq1,lq1);
			// h2=H(ik,jk,lk,inq2,jnq2,lnq2);
			// f+=1;
			// h1+=1;
			// h2+=1;

			f=kernel(iq1,jq1,lq1, inq2, jnq2, lnq2);
			h1=f;
			h2=f;
			f+=1;

			d1re=(float)ft[idq1].re;
			d2re=(float)ft[idq2].re;
			d1im=(float)ft[idq1].im;
			d2im=(float)ft[idq2].im;

			A=d2re*d1re-d2im*d1im;
			//B=d1re*d2re-d1im*d2im;
			//C=d1re*d2re-d1im*d2im;

			//v2[idk].re+=(float)(5.*(h1*A+h2*B)+4.*f*C)/14.;
			v2[idk].re+=(float)(5.*(h1*A+h2*A)+4.*f*A)/14.;

			B=d2re*d1im+d1re*d2im;
			//B=d1re*d2im+d1im*d2re;
			//C=d1re*d2im+d1im*d2re;

			v2[idk].im+=(float)(5.*(h1*B+h2*B)+4.*f*B)/14.;	
			//v2[idk].im+=(float)(5.*(h1*A+h2*B)+4.*f*C)/14.;	

		}

					//t1 = second();

					//cout<< "Time for 1 loop"<< t0-t1<< endl;

	}
	
	
	//t1=MPI_Wtime();
	//t1 = time(NULL);
	t1 = clock();

	cerr<<"loop done"<<endl;

	cout<< "Time: "<< t1-t0<< " " << (t1-t0)/float(idstopn-idstart)<< endl;


	if(part==0){//symmetrize missing values only if all parts are calculated


		for(index=0;index<idstop;index++){

			ik=karr[index*4];
			jk=karr[index*4+1];
			lk=karr[index*4+2];
			idk=karr[index*4+3];

			i=-ik; j=-jk; l=-lk;
			if(i<0){ii=res+i;}else{ii=i;}
			if(j<0){jj=res+j;}else{jj=j;}
			if(l<0){ll=res+l;}else{ll=l;}
			idknew=(ii*res+jj)*res+ll;

	if(idk==0||idk==res/2||idk==res/2*res||idk==res*res/2+res/2||idk==res/2*res*res||idk==res/2*res*res+res/2||idk==(res/2*res+res/2)*res||idk==(res/2*res+res/2)*res+res/2){//cerr<<ik<<" "<<jk<<" "<<lk<<endl; 
	v2[idknew].im=0;}
	else{//continue;
		v2[idknew].re=v2[idk].re;
		v2[idknew].im=-v2[idk].im;}
	}

	output=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_real.txt";
	output2=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_kspace.txt";
	outps=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_ps.txt";

	//smooth with exp(-(kR)^2)
	fftw_complex *ftsm=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	ftsm=Smooth(output.c_str(), output2.c_str(), res, v2, R, Boxsize);

	//calculate power spectrum and output
	Ps(res, Boxsize, 100, ftsm, outps.c_str());

	cerr<<"delta2(x), delta2(k) and ps done"<<endl;

 }


	else{//just output partial calculation

		output2=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_part"+argv[6]+"_kspace.txt";

		ofstream out2(output2.c_str());

		for(i=0;i<res*res*res;i++){out2<<v2[i].re<<" "<<v2[i].im<<endl;}

			out2.close();		

		cerr<<"part " <<argv[6]<< " done"<<endl;

	}
	

	fftw_free(v2);
	free(qarr);
	free(karr);


	return 0;
}
