#include "kernels.hpp" //includes all double/float decisions and FFTW etc
#include "stats.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream> //for cerr
#include <fstream>

int main(int argc, char *argv[]){
	
	if(argc!=7){std::cerr<< "Usage: ./delta2part inputfile res <smoothing scale> <'IC' or 'z0'> Boxsize <part no (0: all, or 1-8)>" <<std::endl; return -1;}

	const std::string arg0=argv[0];
	const std::string argv2=argv[2];
	std::string output, output2, outps;

	const size_t res=atoi(argv[2]);
	const int r2=res/2; //must be int because it will later be compared to ints
	const size_t part=atoi(argv[6]); if(part>8 || part<0){std::cerr<<"part no. not between 0 and 8!"<<std::endl; return -1;}
	size_t idstop=res*res*(r2+1), idstart=(part-1)*(idstop)/8, idstopn=idstart+(idstop)/8; if(part==0){idstart=0; idstopn=idstop;}


	int ik,jk,lk,iq1,jq1,lq1;//,iq2,jq2,lq2;
	size_t idk,idq1,idq2,id,idknew, i, index=0;
	int iq2,jq2,lq2,iq,j,l,ii,jj,ll,iiq2,jjq2,llq2;

	int *qarr=(int*)calloc(res*res*res*4,sizeof(int));
	//int *karr=(int*)calloc(idstop*4,sizeof(int));
	int *karr=(int*)calloc(res*res*res*4,sizeof(int));

	const MyFloat R=atof(argv[3]);	
	const MyFloat Boxsize=atof(argv[5]);
	
	MyFloat d1re, d2re, d1im, d2im,f=0.,h1=0.,h2=0.,A=0.,B=0.;

	fftw_complex *v2=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));

	// //read in delta1(x) from gridfile and ft to delta(k)
	// const char *FILE1=argv[1];
	// size_t RANK=3;

	// hid_t fid1; 
	// hsize_t dims[] = {res,res,res};
	// hid_t datasetid;      // file and memory dataspace identifiers 
	// //hsize_t  offset[RANK];
	
	// dims[0] = res;
	// dims[1] = res;
	// dims[2] = res;

	// //offset[0] = 0;
	// //offset[1] = 0;  
	// //offset[2] = 0;

	// fid1 = H5Fopen(FILE1, H5F_ACC_RDONLY, H5P_DEFAULT);
	// datasetid=H5Dopen(fid1, "Density");
	// hid_t memspacebuf = H5Screate_simple(RANK, dims, NULL );
	// hid_t filespacebuf=H5Dget_space(datasetid);
	// H5Dread(datasetid, H5T_NATIVE_FLOAT, memspacebuf,filespacebuf,H5P_DEFAULT, arr);
	
	// //close identifiers
	// H5Sclose(filespacebuf);
	// H5Sclose(memspacebuf);
	// H5Dclose(datasetid);
	// H5Fclose(fid1);

	const std::string d2file=argv[1];

	MyFloat *in=(MyFloat*)calloc(res*res*res,sizeof(MyFloat));
	std::ifstream d2str(d2file.c_str());
    if (d2str.fail()) {
    	std::cerr << "unable to open file "<<d2file.c_str()<< " for reading" << std::endl;
        exit(1);
    }
    for(i=0;i<res*res*res;i++){d2str>>in[i];}
    d2str.close();

	//fftw_complex *ft = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
    //for(i=0;i<res*res*res;i++){ft[i].re=in[2*i]; ft[i].im=in[2*i+1];}

	//normalise for fftw
	fftw_complex *arr2 = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	for(i=0;i<res*res*res;i++){arr2[i].re=in[i]/pow(res,3.0);}// std::cout<<in[i]*in[i]<<std::endl; } ///pow(res,3.0);}//std::cout<<arr[i]*arr[i]<<std::endl;} 

	//calculate fftw
	fftw_complex *ft = (fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
	fftwnd_plan p = fftw3d_create_plan( res, res, res, FFTW_FORWARD, FFTW_ESTIMATE );
	fftwnd_one(p, arr2, ft);
	fftwnd_destroy_plan(p);
	ft[0].re=0.;
	ft[0].im=0.;

	fftw_free(arr2);
    free(in);

    //for(i=0;i<res*res*res;i++){std::cout << i << " "<< ft[i].re << " " << ft[i].im<<std::endl;}

	//int ll;

	//calulate k values
	for(iq=-(r2-1);iq<r2+1;iq++){
		if(iq<0){ii=res+iq;}else{ii=iq;}
		for(j=-(r2-1);j<r2+1;j++){
			if(j<0){jj=res+j;}else{jj=j;}		
			//for(l=0;l<r2+1;l++){
			for(l=-(r2-1);l<r2+1;l++){
				if(l<0){ll=res+l;}else{ll=l;}
				
				idk=(ii*(res)+jj)*(res)+ll;		
				
				karr[4*index]=iq; karr[4*index+1]=j; karr[4*index+2]=l; karr[4*index+3]=idk;
				index+=1;

			}
		}
	}


	index=0;
	//calculate q values
	for(iq=-(r2-1);iq<r2+1;iq++){
		if(iq<0){ii=res+iq;}else{ii=iq;}
		for(j=-(r2-1);j<r2+1;j++){
			if(j<0){jj=res+j;}else{jj=j;}		
			//for(l=-(r2-1);l<r2+1;l++){
			for(l=0;l<r2+1;l++){	
				if(l<0){ll=res+l;}else{ll=l;}
				
				idk=(ii*(res)+jj)*(res)+ll;				
				qarr[4*index]=iq; qarr[4*index+1]=j; qarr[4*index+2]=l;  qarr[4*index+3]=idk;
				index+=1;
			}
		}
	}


	std::cerr<<"beginning loop"<<std::endl;

	int lowi, lowj, lowl, hii, hij, hil;
	int iiq1, jjq1, llq1, iqn2, jqn2, lqn2;

	clock_t t0, t1;
	t0=clock();

	idstopn=res*res*res;
	//part=15;

	for(id=idstart;id<idstopn;id++){

		ik=karr[id*4];
		jk=karr[id*4+1];
		lk=karr[4*id+2];
		idk=karr[4*id+3];

		//std::cout << id << " " << idk << std::endl;

		lowi=-r2+1;
		lowj=-r2+1;
		lowl=-r2+1;		

		hii=r2+1;
		hij=r2+1;
		hil=r2+1;

		// lowi=std::max(ik-r2, -r2);
		// lowj=std::max(jk-r2, -r2);
		// lowl=std::max(lk-r2, -r2);		

		// hii=std::min(ik+r2, r2)+1;
		// hij=std::min(jk+r2, r2)+1;
		// hil=std::min(lk+r2, r2)+1;

		if(idk == 22 || idk ==62){std:: cout << idk<< ": "<< lowi << " " << hii << " " << lowj << " " << hij << " " << lowl << " " << hil << std::endl;}

		for(iq1 = lowi; iq1<hii; iq1++){
			if(iq1<0){iiq1=res+iq1;}else{iiq1=iq1;}
			iq2=ik-iq1;
			if(iq2<0){iiq2=res+iq2;}else{iiq2=iq2;} //never goes > res-1, so two conditions are enough
			
			for(jq1 = lowj; jq1<hij; jq1++){
				if(jq1<0){jjq1=res+jq1;}else{jjq1=jq1;}
				jq2=jk-jq1;
				if(jq2<0){jjq2=res+jq2;}else{jjq2=jq2;}

				for(lq1 = lowl; lq1<hil; lq1++){
					if(lq1<0){llq1=res+lq1;}else{llq1=lq1;}
					lq2=lk-lq1;
					if(lq2<0){llq2=res+lq2;}else{llq2=lq2;}									

					idq1=(iiq1*res+jjq1)*res+llq1;				

					if(iq2<-r2+1){iqn2=res+iq2;}else if(iq2>r2+1){iqn2=iq2-res;}else{iqn2=iq2;} //never goes > res-1, so two conditions are enough
					if(jq2<-r2+1){jqn2=res+jq2;}else if(jq2>r2+1){jqn2=jq2-res;}else{jqn2=jq2;}
					if(lq2<-r2+1){lqn2=res+lq2;}else if(lq2>r2+1){lqn2=lq2-res;}else{lqn2=lq2;}		

					idq2=(iiq2*res+jjq2)*res+llq2;

					//kernels:
					// f=beta(iq1,jq1,lq1,iq2,jq2,lq2);
					// h1=alpha(ik,jk,lk,iq1,jq1,lq1);
					// h2=alpha(ik,jk,lk,iq2,jq2,lq2);
					
					 f=beta(iq1,jq1,lq1,iqn2,jqn2,lqn2);
					 h1=alpha(ik,jk,lk,iq1,jq1,lq1);
					 h2=alpha(ik,jk,lk,iqn2,jqn2,lqn2);
					
					//f=1.;
					//h1=1.;
					//h2=1.;

					//f=qarr[idk];
					//h1=qarr[idq1];
					//h2=qarr[idq2];

					d1re=(MyFloat)ft[idq1].re;
					d2re=(MyFloat)ft[idq2].re;
					d1im=(MyFloat)ft[idq1].im;
					d2im=(MyFloat)ft[idq2].im;

					A=d2re*d1re-d2im*d1im;
					v2[idk].re+=(5.*(h1+h2)+4.*f)/14.*A;

					B=d2re*d1im+d1re*d2im;
					v2[idk].im+=(5.*(h1+h2)+4.*f)/14.*B;		

					// if(idk == 22 || idk ==62){ 
					//  	std::cout << idk << " (" << ik << ", "<< jk << ", "<< lk << ") | "<< h1 << " " << h2 << " "<< f <<  " = " << (5*h1+5*h2+4*f)/14. << std::endl;
					//  	std::cout << idq1 << "( "<<ft[idq1].re<< ", "<<ft[idq1].im<< " ); " << idq2 <<"( "<<ft[idq2].re<< ", "<<ft[idq2].im<< " ) " << std::endl;
					//  	std::cout << "("<< iq1 << ", " << jq1 << ", " << lq1 << "); " <<  "("<< inq2 << ", " << jnq2 << ", " << lnq2 << ")" << "| A: " << (5.*(h1+h2)+4.*f)/14.*A << " |  B: " <<(5.*(h1+h2)+4.*f)/14.*B << std::endl;
					// // 	std::cout << h1 << " "<< h2 << " " << f << ", "<< A << " " << B <<std::endl;
					// // 	std::cout << std::endl;
					// }


			
				}
			}
		}
		
		//t1 = second();
		//std::cout<< "Time for 1 loop"<< t0-t1<< std::endl;
	
	} //end of idk loop
	
	t1 = clock();	

	std::cerr<<"loop done"<<std::endl;
	std::cout<< "Time: "<< t1-t0<< " " << (t1-t0)/MyFloat(idstopn-idstart)<< std::endl;


	if(part==0){//symmetrize missing values only if all parts are calculated
		 v2[r2].im=0.;
		 v2[r2*res].im=0;
		 v2[r2*res+r2].im=0;
		 v2[r2*res*res].im=0;
		 v2[r2*res*res+r2].im=0;
		 v2[(r2*res+r2)*res].im=0;
		 v2[(r2*res+r2)*res+r2].im=0;
		int index2=0;
		//std:: cout << (res*res*res)-idstopn << std::endl;
		for(index=0;index<idstop;index++){ //index=0 is dealt with below
		//for(index=0;index<(res*res*res)-idstopn;index++){
		//for(index=0;index<res*res*res;index++){
		// for(iq=-(r2-1);iq<r2+1;iq++){
		// 	if(iq<0){ii=res+iq;}else{ii=iq;}
		// 	for(j=-(r2-1);j<r2+1;j++){
		// 		if(j<0){jj=res+j;}else{jj=j;}		
		// 		for(l=-(r2-1);l<0;l++){
		// 		if(l<0){ll=res+l;}else{ll=l;}


			//index2=(index*(res-1)+1)*4; //was correct for res=4 at some point but now it's wrong
			//std::cout << index << " " << index2 << " " << std::endl;
			// ik=karr[index2];
			// jk=karr[index2+1];
			// lk=karr[index2+2];
			// idk=karr[index2+3];

			ik=qarr[index*4];
			jk=qarr[index*4+1];
			lk=qarr[index*4+2];
			idk=qarr[index*4+3];
			//ik=-iq; jk=-j; lk=-l;

			iq=-ik; j=-jk; l=-lk;
			if(iq<0){ii=res+iq;}else{ii=iq;}
			if(j<0){jj=res+j;}else{jj=j;}
			if(l<0){ll=res+l;}else{ll=l;}
			idknew=(ii*res+jj)*res+ll;

			//if( fabs(v2[idk].im + v2[idknew].im ) >1e-7 ){
			//std::cout << index << " "<< idk << " " << idknew << " (" << ik << ", " << jk << ", "<< lk << "); (" << iq << ", " << j << ", "<< l<< ") | " <<  v2[idk].im << " " << v2[idknew].im << " "<< fabs(v2[idk].im + v2[idknew].im )<< " | " << v2[idk].re << " " << v2[idknew].re <<std::endl;
			//std:: cout << index << " "<< idk << " " << idknew << " (" << ik << ", " << jk << ", "<< lk << "); (" << iq << ", " << j << ", "<< l<< ")" << std::endl;
			//std::cout << std::endl;
			//}
			//std:: cout << index << " "<< idk << " " << idknew << " (" << v2[idk].re << ", " << v2[idk].im << "); (" << v2[idknew].re << ", " << v2[idknew].im<< ")" << std::endl;
			//std::cout << std::endl;

			v2[idknew].re=v2[idk].re;
			v2[idknew].im=-v2[idk].im;
		
		}

		v2[0].re=0.; //enforce the mean 0 condition (already fine up to numerical accuracy anyway)
		v2[0].im=0.;

		output=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_real.txt";
		output2=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_kspace.txt";
		outps=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_ps.txt";

		//smooth with exp(-(kR)^2)
		fftw_complex *ftsm=(fftw_complex*)calloc(res*res*res,sizeof(fftw_complex));
		ftsm=Smooth(output.c_str(), output2.c_str(), res, v2, R, Boxsize);

		//calculate power spectrum and output
		Ps(res, Boxsize, 100, ftsm, outps.c_str());

		std::cerr<<"delta2(x), delta2(k) and ps done"<<std::endl;

 	}


	else{//just output partial calculation

		output2=arg0+"_"+argv[2]+"_"+argv[4]+argv[5]+"_"+argv[3]+"_part"+argv[6]+"_kspace.txt";

		std::ofstream out2(output2.c_str());
		for(i=0;i<res*res*res;i++){out2<<v2[i].re<<" "<<v2[i].im<<std::endl;}
		out2.close();		

		std::cerr<<"part " <<argv[6]<< " done"<<std::endl;

	}
	

	fftw_free(v2);
	free(qarr);
	free(karr);

	return 0;
}
