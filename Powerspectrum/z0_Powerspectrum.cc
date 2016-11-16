#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <error.h>

#include <srfftw.h>
#include "HDF_IO.hh"

bool DoesFileExist( char *Filename ){
        bool flag = false;
        std::fstream fin(Filename,std::ios::in);
        if( fin.is_open() )
                flag=true;
        fin.close();
        return flag;
}

void SaveHDF( char Filename[], int m_nGrid[3], float Boxlength, float *m_Data )
{
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t       offset[3], count[3];

  
  file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  for( int i=0; i<3; ++i )
        count[i] = m_nGrid[2-i];

  filespace = H5Screate_simple( 3, count, NULL );
  dset_id = H5Dcreate( file_id, "Density", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT );
  H5Sclose(filespace);

  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;

 // std::cout<< "5234 "<< data2[0]<<" "<< data2[1]<<" "<<data2[2]<< std::endl;
 // std::cout<< "786 " <<m_Data[0]<<" "<< m_Data[1]<<" "<<m_Data[2]<< std::endl;

  memspace = H5Screate_simple( 3, count, NULL );

  filespace= H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, m_Data);
  H5Sclose(memspace);
  H5Dclose(dset_id);
  H5Sclose(filespace);


   H5Fclose(file_id);
}

void InitGrid( int nx, int ny, int nz, fftw_real *data )
{
unsigned ix, iy, iz;
  
  for( int ix=0; ix<nx; ++ix )
    for( int iy=0; iy<ny; ++iy )
      for( int iz=0; iz<nz; ++iz )
        data[(ix*ny + iy) * nz + iz] = -1.0f;
}

void GridParticles( int nx, int ny, int nz, std::vector<float> &Pos, float Boxlength, unsigned nPart, float m,  fftw_real *data , float subbox, int *ourbox)
{
  
    unsigned ix, iy, iz, i;

  //double wpar = 1.0f/((double)nx*ny*nz)/(double)nPart;//pow(Boxlength,3)/(double)nPart;//((double)nx*ny*nz)/(double)nPart;
  float wpar =  (float)(nx*ny*nz)/(float)nPart;
  /*std::cerr << "lbox = " << Boxlength << std::endl;
  std::cerr << "nPart= " << nPart << std::endl;
  std::cerr << "wpar = " << wpar << std::endl;*/
  float x,y,z,dx,dy,dz,tx,ty,tz,tyw,dyw;
  unsigned ix1,iy1,iz1;
  int icount = 0;
  unsigned nPartThis = Pos.size()/3;
  //int ourbox[3]  = {4, 4, 4};

//  printf("test ourbox = %d %d %d \n", ourbox[0], ourbox[1], ourbox[2]);


  for (i = 0; i < 3; i++)
      if (ourbox[i] > subbox || ourbox[i] < 1)
	  error (-1,0,"Error with ourbox.");

  // bad approx    wpar *= subbox;   // tag+
     

  for( unsigned i=0; i<nPartThis; ++i ){
    x = Pos[3*i+0];
    y = Pos[3*i+1];
    z = Pos[3*i+2];

// if x,y,z < nx,ny,nz :

    x -= (ourbox[0]-1.) * (float)nx;
    y -= (ourbox[1]-1.) * (float)ny;
    z -= (ourbox[2]-1.) * (float)nz;

      if (x < nx && y < ny && z < nz && x > 0. && y > 0. && z > 0.){
	  ++icount;
      }
  }

//  printf("icount = %d \n", icount);

  wpar =  (float)(nx*ny*nz)/(float)nPart;

// (float)(nx*ny*nz)/(float)icount;

  for( unsigned i=0; i<nPartThis; ++i ){
  
    x = Pos[3*i+0];
    y = Pos[3*i+1];
    z = Pos[3*i+2];
  

// if x,y,z < nx,ny,nz :

    x -= (ourbox[0]-1.) * (float)nx;
    y -= (ourbox[1]-1.) * (float)ny;
    z -= (ourbox[2]-1.) * (float)nz;

    //printf("x,y,z = %f \t %f \t %f\n", x, y, z);

    if (x < nx && y < ny && z < nz && x > 0. && y > 0. && z > 0.){

//	++icount;

    ix = (unsigned)x;
    iy = (unsigned)y;
    iz = (unsigned)z;

    dx = (x-((float)ix));
    dy = (y-((float)iy));
    dz = (z-((float)iz));

    ix %= nx;
    iy %= ny;
    iz %= nz;

    tx = 1.0f-dx;
    ty = 1.0f-dy;
    tz = 1.0f-dz;

    tyw = ty*wpar;
    dyw = dy*wpar;

    ix1 = (ix+1)%nx;
    iy1 = (iy+1)%ny;
    iz1 = (iz+1)%nz;
    
    data[(ix*ny + iy) * nz + iz]   += tz*tx*tyw;
    data[(ix1*ny + iy) * nz + iz]  += tz*dx*tyw;
    data[(ix*ny + iy1) * nz + iz]  += tz*tx*dyw;
    data[(ix1*ny + iy1) * nz + iz] += tz*dx*dyw;

    data[(ix*ny + iy) * nz + iz1]   += dz*tx*tyw;
    data[(ix1*ny + iy) * nz + iz1]  += dz*dx*tyw;
    data[(ix*ny + iy1) * nz + iz1]  += dz*tx*dyw;
    data[(ix1*ny + iy1) * nz + iz1] += dz*dx*dyw;
    }

  }
}

//void ComputePowerSpectrum( char Filename[], int nBins, int nx, int ny, int nz, float Boxlength, int nPart, fftw_real *data )
void ComputePowerSpectrum( char Filename[], int nBins, int nx, int ny, int nz, float Boxlength, fftw_real *data )
{
  fftw_complex	*in;
  fftwnd_plan 	p;
  int		ix, iy, iz, idx, idx2,idx3;
  unsigned   *inBin = new unsigned[nBins];
  unsigned   *inBin2 = new unsigned[nBins];
  float       k, kfft,
    //kmax = M_PI/Boxlength*nx,
    kmax = M_PI,
    kmin = 2.0f*M_PI/(float)nx,
    dklog = log10(kmax/kmin)/nBins,
   // dk=(0.4*Boxlength/(float)nx-kmin)/nBins,
    kw = 2.0f*M_PI/(float)nx;

	//std::cout<<dk<<std::endl;
    
  // test break box tag+
//  nx /= 2.;  // tag+
//  ny /= 2.;  // tag+
//  nz /= 2.;  // tag+
  
  //fftw_complex ft[nx][ny][(nz/2+1)];
  fftw_complex *ft = new fftw_complex[nx*ny*(nz/2+1)];
  p = rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );//|FFTW_IN_PLACE );
  rfftwnd_one_real_to_complex(p, data, ft );
  in = (fftw_complex*) data;
  
  float *Gx = new float[nBins];
  float *kbin = new float[nBins];
 float *Gx2 = new float[nBins];
  float *kbin2 = new float[nBins];
  for( ix=0; ix<nBins; ix++ ){
    inBin[ix] = 0;
    Gx[ix] = 0.0f;
    kbin[ix] = 0.0f;
    inBin2[ix] = 0;
    Gx2[ix] = 0.0f;
    kbin2[ix] = 0.0f;
  }
  
	std::ofstream file;
	file.open("P_pdf.txt");  

  for(ix=0; ix<nx;ix++)
    for(iy=0;iy<ny;iy++)
      for(iz=0;iz<nz/2+1;iz++){
        idx = (ix*ny+iy)*(nz/2+1)+iz;

        // determine mode modulus
        //float vabs = ft[ix][iy][iz].re*ft[ix][iy][iz].re+ft[ix][iy][iz].im*ft[ix][iy][iz].im;
	float vabs = ft[idx].re*ft[idx].re+ft[idx].im*ft[idx].im;
        int iix, iiy, iiz;

        if( ix>nx/2 ) iix = ix - nx; else iix = ix;
        if( iy>ny/2 ) iiy = iy - ny; else iiy = iy;
        iiz = iz;
	
        kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);
        k = kfft*kw;
	
        // correct for aliasing, formula from Jing (2005), ApJ 620, 559
        // assume isotropic aliasing (approx. true for k<kmax=knyquist)
        // this formula is for CIC interpolation scheme, which we use
        float JingCorr = (1.0f-2.0f/3.0f*sin(M_PI*k/kmax/2.0f)*sin(M_PI*k/kmax/2.0f));	
	vabs /= JingCorr;
        	
        //.. logarithmic spacing in k
	      idx2 = (int)(( 1.0f/dklog * log10(k/kmin) ));
	      //idx3 = (int)((k-kmin)/dk);	

        if(k>=kmin&&k<kmax){
	if(idx2==15)
		file << ix << "\t" << iy << "\t" << iz << "\t" << vabs << std::endl;
 
        
	  if( iz == 0 ){
            Gx[idx2] += vabs;
            kbin[idx2] += k;
            inBin[idx2]++;
// 	    Gx2[idx3] += vabs;
//             kbin2[idx3] += k;
//             inBin2[idx3]++;
          }else {
            Gx[idx2] += 2.0f*vabs;
            kbin[idx2] += 2.0f*k;
            inBin[idx2]+= 2;
// 	    Gx2[idx3] += 2.0f*vabs;
//             kbin2[idx3] += 2.0f*k;
//             inBin2[idx3]+= 2;
          }
        }
      }
    

	file.close();
  //... convert to physical units ...
  std::ofstream ofs(Filename);
//std::ofstream ofs2("test.dat");
  
  // FFT norm brings 1/N^3 (actually it gives delta/N^(3/2), but Powersp=delta*delta), conversion to physical units L^3/N^3, definition of powerspectrum brings (2pi)^-3
  float fftnorm = pow(Boxlength,3.0)/pow(2.0*M_PI,3.0);
  
  for(ix=0; ix<nBins; ix++){
	    
	if(inBin[ix]>0){
    ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) + log10 ((float)nx/(float)Boxlength))
	<< std::setw(16) <<  kbin[ix]/inBin[ix]*(float)nx/(float)Boxlength
        << std::setw(16) << ((double)((Gx[ix]/inBin[ix])/nx/ny/nz/nx/ny/nz))*fftnorm
        << std::setw(16) << inBin[ix]
        << std::endl;
  }
// 	else if(inBin2[ix]>0){ofs2 <<std::setw(16) << (kmin+dk*(ix+0.5))/(float)Boxlength*(float)nx
// 	<< std::setw(16) << ((double)(Gx2[ix]/inBin2[ix]/nx/ny/nz/nx/ny/nz))*fftnorm
// 	<< std::setw(16) << inBin2[ix]<< std::endl;
// 	}

	else{continue;}
   }
  
  fftwnd_destroy_plan(p);
  delete[] inBin;
  delete[] Gx;
  delete[] kbin;

  delete[] ft;
}


//computation of the bispectrum in a very simple configuration: K, K, -2k.

void ComputeBispectrum( char Filename[], int nBins, int nx, int ny, int nz, float Boxlength, int nPart, fftw_real *data )
{
  fftw_complex	*in;
  fftwnd_plan 	p;
  int		ix, iy, iz, idx, idx2;
  int 		tidx;
  unsigned   *inBin = new unsigned[nBins];
  float      
    k, kfft,
    //kmax = M_PI/Boxlength*nx,
    kmax=M_PI,
    kmin = 2.0f*M_PI/(float)nx,
    dklog = log10(kmax/kmin)/nBins,
    kw = 2.0f*M_PI/(float)nx;
      
  
  //fftw_complex ft[nx][ny][(nz/2+1)];
  fftw_complex *ft = new fftw_complex[nx*ny*(nz/2+1)];
  p = rfftw3d_create_plan( nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );//|FFTW_IN_PLACE );
  rfftwnd_one_real_to_complex(p, data, ft );
  in = (fftw_complex*) data;
  
  double *Gx = new double[nBins];
  double *sGx = new double[nBins]; //sigma of the Bispectrum	
  float *kbin = new float[nBins];

  for( ix=0; ix<nBins; ix++ ){
    inBin[ix] = 0;
    Gx[ix] = 0.0;
    sGx[ix] = 0.0;
    kbin[ix] = 0.0f;
  }

	std::ofstream file;
	file.open("pdf.txt");	  
  
  for(ix=-nx/4+1; ix<nx/4;ix++)
    for(iy=-ny/4+1;iy<ny/4;iy++)
      //for(iz=0;iz<nz/2+1;iz++)
	for(iz=0;iz<nz/4;iz++)	
	{
	int iix,iiy;

	iix = (ix+nx)%nx;
	iiy = (iy+ny)%ny;
        idx = (iix*ny+iiy)*(nz/2+1)+iz;
	iix = (2*ix+nx)%nx;
	iiy = (2*iy+ny)%ny;
 	tidx= (iix*ny+iiy)*(nz/2+1)+2*iz;

        kfft = sqrt(ix*ix+iy*iy+iz*iz);
        k= kfft*kw;
	
	//std::cout << ix << "\t" << iy << "\t" << k << std::endl;
	float reValue = ft[idx].re*ft[idx].re*ft[tidx].re-ft[idx].im*ft[idx].im*ft[tidx].re+2.*ft[idx].re*ft[idx].im*ft[tidx].im;
	float squaredReValue = powf(reValue,2.);
	//float imValue = 2.*ft[idx].re*ft[idx].im*ft[tidx].re-ft[idx].re*ft[tidx].im+ft[idx].im*ft[idx].im*ft[tidx].im;
        //the imaginary part is necessarily zero, taking into account also the fourth of the cube which 
        //.. logarithmic spacing in k
	      idx2 = (int)(( 1.0f/dklog * log10(k/kmin) ));
	
	//std::cout << "kmax= " << kmax << " kmin= " << kmin << " k= " <<k<< std::endl;
        if(k>=kmin&&k<kmax){
	//std::cout << "kmax= " << kmax << " kmin= " << kmin << " k= " <<k<< std::endl;
	//std::cout << idx2<< std::endl;
		if(idx2==15)
		file << reValue << std::endl;


  
  // FFT norm brings 1/N^3, conversion to physical units (L/N/2pi)^6, definition of bispectrum brings (2pi)^-6
  float fftnorm = pow(Boxlength,3.0)/pow(2.0*M_PI,3.0);

	  if( iz == 0 ){
            Gx[idx2] += ((double)(reValue)/nx/ny/nz/nx/ny/nz/nx/ny/nz*pow(fftnorm,2.));
	    sGx[idx2] += (double)(squaredReValue)/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz*powf(fftnorm,4.);
            kbin[idx2] += k;
            inBin[idx2]++;
          }else {
            Gx[idx2] += 2.0*(double)(reValue)/nx/ny/nz/nx/ny/nz/nx/ny/nz*pow(fftnorm,2.);
	    sGx[idx2] += 2.0*(double)(squaredReValue)/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz/nx/ny/nz*powf(fftnorm,4.);
            kbin[idx2] += 2.0f*k;
            inBin[idx2]+= 2;
          }
        }
	}	
	file.close();


  std::ofstream ofs(Filename);

  
  for(ix=0; ix<nBins; ix++){
    
      ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) + log10 ((float)nx/(float)Boxlength))
	<< std::setw(16) <<   kbin[ix]/inBin[ix]*(float)nx/(float)Boxlength
        << std::setw(16) << ((double)((Gx[ix]/inBin[ix])))
	<< std::setw(16) << powf(((double)((sGx[ix]/(inBin[ix]-1))))-((double)(powf((Gx[ix]/inBin[ix]),2.))),0.5)
        << std::setw(16) << inBin[ix]
        << std::endl;
  }

  
  fftwnd_destroy_plan(p);
  delete[] inBin;
  delete[] Gx;
  delete[] sGx;
  delete[] kbin;

  delete[] ft;
}







int main( int argc, char *argv[])
{

  std::vector<float> Pos;

  if( argc != 8 ){
    std::cerr << " Usage: Powerspectrum <gadgetdata.hdf5> <ngrid> <spectrumname.dat> <gridfilename.dat> ourbox x,y,z \n" << std::endl;
    return -1;
  }

  float subbox = 1.; // tag+

  float Boxlength;
  try{
    HDFReadGroupAttribute( argv[1], "/Header", "BoxSize", Boxlength );
  }catch(...){
    std::cout << " Please specify box size in Mpc/h: ";
    std::cin >> Boxlength;  
  }

  Boxlength /= subbox; // tag+

  std::cout << " Box is " << Boxlength << std::endl;
  
  std::vector<float> masses;
  HDFReadGroupAttribute( argv[1], "/Header", "MassTable", masses );
  float pmass(masses[1]);
 
  //  float pmass;
  //  std::cout << " Please specify particle mass in 1e10 M_sun/h: ";
  //  std::cin >> pmass;  
  std::cout << " Particle mass is " << pmass << std::endl;

  std::vector<unsigned> numpart;
  HDFReadGroupAttribute( argv[1], "/Header", "NumPart_Total", numpart );
  unsigned numpartDM(numpart[1]);
  //   unsigned numpartDM;	  
  //   std::cout << " Please specify total number of particle: ";
  //   std::cin >> numpartDM;  
   
     std::cout << " found " << numpartDM << " particles\n";

  int ngrid = atoi(argv[2]);
  //std::cout << " P											lease specify grid resolution for FFT: ";
  //std::cin >> ngrid;
  
  int ourbox[3];
  ourbox[0] = atoi(argv[5]); //atoi: converts (input)string to integer number (boxsize in mpc/h)
  ourbox[1] = atoi(argv[6]);
  ourbox[2] = atoi(argv[7]);

  printf("ourbox = %d %d %d \n", ourbox[0], ourbox[1], ourbox[2]);
 
  int nx=ngrid, ny=ngrid, nz=ngrid;

  fftw_real *data = new fftw_real[nx*ny*nz];

  InitGrid(nx,ny,nz,data);
//   std::cout<< "4123 " <<data[0]<<" "<<data[1]<< " "<< data[2]<< std::endl;

  unsigned nperslab = 10000000, nread=0;

  
  while(nread < numpartDM){
	if( nread+nperslab>numpartDM )
		nperslab = numpartDM-nread;
  	HDFReadVectorSlab( argv[1], "/PartType1/Coordinates", nread, nperslab, Pos );
	nread += nperslab;
//	std::cerr << "Read " << nread << " of " << numpartDM << ".\n";
  	for( unsigned i=0; i<Pos.size(); ++i ){
    		Pos[i] *= (float)ngrid/Boxlength;
  	}
  	//unsigned npart = Pos.size()/3;

  	GridParticles( nx, ny, nz, Pos, Boxlength, numpartDM,  pmass, data, subbox, ourbox );
  }

  int nres[3] = {nx, ny, nz};
 // std::cout<< "1234 "<< data[0]<<" "<<data[1]<< " "<< data[2]<< std::endl;
  SaveHDF( argv[4], nres, Boxlength, (float*)data );

  int nBins = 100;//15;//(int)(ngrid/1.0);
 
  std::cout << "nBins= " << nBins << std::endl;
  //ComputePowerSpectrum( argv[3], nBins, nx, ny, nz, Boxlength, numpartDM, data );
  ComputePowerSpectrum( argv[3], nBins, nx, ny, nz, Boxlength, data );

  //ComputeBispectrum(argv[4], nBins, nx, ny, nz, Boxlength, numpartDM, data );
  
  delete[] data;
  return 1;
}

