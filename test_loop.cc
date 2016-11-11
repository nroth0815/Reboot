#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream> //for cerr
#include <fstream>

#include <rfftw.h>

int main(int argc, char *argv[]){
	
	if(argc!=2){std::cerr<< "Usage: ./delta2part res " <<std::endl; return -1;}

	//const size_t res=atoi(argv[1]);
	int iq, jq, lq, index=0;
	//clock_t t0, t1;
	int i;
	int res=8;
	
	int *qarr=(int*)calloc(res*res*res*4,sizeof(int));
	int q;

	double sinarr[8]={ 1., 7.07106781e-01, 0, -7.07106781e-01,  -1,  -7.07106781e-01, 0,   7.07106781e-01};
	//{ 0,   7.07106781e-01,   1., 7.07106781e-01, 0,  -7.07106781e-01, -1.,  -7.07106781e-01}; //sine
	//{ 1., 7.07106781e-01, 0, -7.07106781e-01,  -1,  -7.07106781e-01, 0,   7.07106781e-01}; //cosine
	

    for(iq=0; iq< res ; iq++){
    	std::cout << "sinarr: "<< sinarr[iq] << std::endl;
    }

    std::cout <<std::endl;

    fftw_complex *arr2 = (fftw_complex*)calloc(res, sizeof(fftw_complex));
	for(i=0;i<res;i++){arr2[i].re=sinarr[i]; } ///pow(res,3.0);}//std::cout<<arr[i]*arr[i]<<std::endl;} 

	//calculate fftw
	fftw_complex *ft = (fftw_complex*)calloc(res, sizeof(fftw_complex));
	fftw_plan p = fftw_create_plan( res, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_one(p, arr2, ft);
	fftw_destroy_plan(p);
	ft[3].re=0.; //get rid of numerical inaccuaries in sampled sine that translate to fourier space
	ft[5].re=0.; //get rid of numerical inaccuaries in sampled sine that translate to fourier space

    for(iq=0; iq< res ; iq++){
    	std::cout << "ft(sinarr): "<< ft[iq].re << " " << ft[iq].im << std::endl;
    }

	fftw_complex *arr3 = (fftw_complex*)calloc(res, sizeof(fftw_complex));

	int iiq, jjq, llq;
	fftw_complex kernel={1.,0};
	//fftw_complex test=ft[0]*kernel2;
	//std:: cout << test.re << " " << test.im<< std::endl;
	//double kernel=1.;

	for(iq=-res/2+1; iq< res/2+1; iq++){
		for(jq=-res/2+1; jq< res/2+1; jq++){
			lq=iq-jq;

			if(iq < 0 ){iiq=iq+res;}
			else{iiq=iq;}

			if(jq < 0 ){jjq=jq+res;}
			else{jjq=jq;}

			if(lq < 0 ){llq=lq+res;}
			else{llq=lq;}

			kernel.re=0.;

			//if(lq < (-res/2+1)){kernel.im = lq+res;}
			//else if(lq > (res/2) ){kernel.im = lq-res;}
			//else{kernel.im=llq;} //wrapped -> gives complex-valued result! -> WRONG
			//kernel.im=iiq; //convolution in first variable, wrapped -> gives complex-valued result! -> WRONG

			kernel.im=lq; //not wrapped; gives correct result!
			//kernel.im=iq; //not wrapped; convolution with first argument; correct up to factor 2!

    		arr3[iiq].re+= (ft[jjq].re * ft[llq].re - ft[jjq].im * ft[llq].im)*kernel.re +
    		 (ft[jjq].re * ft[llq].im + ft[jjq].im * ft[llq].re)*kernel.im;
    		arr3[iiq].im+= (ft[jjq].re * ft[llq].im + ft[jjq].im * ft[llq].re)*kernel.re +
    		 (ft[jjq].re * ft[llq].re - ft[jjq].im * ft[llq].im)*kernel.im;
    	}
    }

	std::cout << std::endl;

    for(iq=0; iq< res ; iq++){
    	std::cout << "summed: "<< arr3[iq].re << " " << arr3[iq].im << std::endl;
    }

    fftw_complex *out = (fftw_complex*)calloc(res, sizeof(fftw_complex));
	fftw_plan p2 = fftw_create_plan( res, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_one(p2, arr3, out);
	fftw_destroy_plan(p2);

	std::cout << std::endl;

	for(iq=0; iq< res ; iq++){
    	std::cout << "output: " << out[iq].re/res/res << " " << out[iq].im/res/res << std::endl;
    }

	// t0=clock();

	// for(iq=0;iq<res;iq++){
	// 	for(jq=0;jq<res;jq++){
	// 		for(lq=0;lq<res;lq++){

	// 			index+=1;
	// 			q=qarr[index];
	// 			q+=1;

	// 		}
	// 	}
	// }


	// t1=clock();

	// std::cerr<<"loop done"<<std::endl;
	// std::cout<< "Total time: "<< t1-t0<< ",  per iteration: " << (t1-t0)/float(res*res*res)<< std::endl;


	free(qarr);
	fftw_free(arr2);
	fftw_free(ft);

	return 0;

}


