#include <string>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <boost/functional/hash.hpp> //better hash libraries (for later)
#include "hash_temp.hpp" //contains Key and Keyhasher*

#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions

#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

#include <limits>



struct ab{

    MyFloat alpha;
    MyFloat beta;

} values;

void calc_ab(size_t res, ab *array){

    for(size_t i=0; i < res*res*res; i++){
        array[i].alpha=2*i;
        array[i].beta=-i;
    }

}

void test_hash(ab *array, size_t res, gsl_rng * r){

    typedef std::unordered_map< Key, ab, KeyHasher> hashmap;
    hashmap numbers;

    Key index={0,0};
    for(size_t i=0; i<(res*res*res); i++){
        index.idq=(size_t)(gsl_ran_flat(r, 0.,1.)*(res*res*res)-1);
        numbers[index]=array[i];
    }

    double test=0;    
    std::chrono::time_point<std::chrono::system_clock> startP,endP;
    startP = std::chrono::system_clock::now();
     //hashmap::hasher hashfunc= numbers.hash_function();
    for(size_t i=0; i<(res*res); i++){
        index.idq=(size_t)(gsl_ran_flat(r, 0.,1.)*(res*res*res)-1);
        //index.idp=(size_t)(gsl_ran_flat(r, 0.,1.)*(res*res*res)-1);
        //std::cout << index.idq << " "<< index.idp << std::endl;
        test+=numbers[index].alpha + numbers[index].beta;

    }

    // for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
    //      test+= i->second.alpha + i->second.beta;
          //std::cout << i->first.idq << ", " << i->first.idp << " -> "<< i->second.alpha << " " << i->second.beta << "; (hash = " << hashfunc( i->first ) << ")"  << std::endl;
     //}
     endP = std::chrono::system_clock::now();
     std::cout << "duration hash-loop (seconds):  " << std::chrono::duration_cast<std::chrono::nanoseconds>(endP-startP).count()/1e09 << std::endl;
     std::cout << test/1e8 << std::endl;
     std::cout << numbers.bucket_count() << " " << numbers.max_bucket_count() << " " << std::endl;
     std::cout << numbers.load_factor() << " " << numbers.max_load_factor() << std::endl;


}

void test_array(ab *array, size_t res, gsl_rng * r){

    ab *ab_arr2=(ab*)calloc(res*res*res,sizeof(ab));

    size_t index;
    for(size_t i=(res*res*res)-1; i>0; i--){
        index=(size_t)(gsl_ran_flat(r, 0.,1.)*(res*res*res)-1);
        ab_arr2[index]=array[index];
    }


    double test=0;

    std::chrono::time_point<std::chrono::system_clock> startP,endP;
    startP = std::chrono::system_clock::now();
    for(size_t i=0; i<(res*res); i++){
        index=(size_t)(gsl_ran_flat(r, 0.,1.)*(res*res*res)-1);
        //std::cout << index << std::endl;
        //index=gsl_ran_gaussian_ziggurat(r,1.);
        //index.idp=gsl_ran_gaussian_ziggurat(r,1.);
        test+=ab_arr2[index].alpha+ab_arr2[index].beta;

    }
    // for(size_t i=(res*res*res)-1; i>0; i--){
    //     test+=ab_arr2[i].alpha + ab_arr2[i].beta;
    // }
    endP = std::chrono::system_clock::now();
    std::cout << "duration array-loop (seconds):  " << std::chrono::duration_cast<std::chrono::nanoseconds>(endP-startP).count()/1e09 << std::endl;
    std::cout << test/1e8<< std::endl;

}

int main(int argc, char *argv[]){

    if(argc!=2){std::cerr << "Only one argument expected!" << std::endl; exit(1);}

    size_t res=atoi(argv[1]);
    
    gsl_rng * r;
    const gsl_rng_type * T;
    int seed =342345110;

#ifndef DOUBLEPRECISION
T = gsl_rng_ranlxs2; //this is the name of the generator defined in the environmental variable GSL_RNG_TYPE
#else
T = gsl_rng_ranlxd2; 
#endif

    r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
    gsl_rng_set(r,seed);

    ab *ab_arr=(ab*)calloc(res*res*res,sizeof(ab));

    calc_ab(res, ab_arr);

    uint64_t duration1, duration2;

    std::chrono::time_point<std::chrono::system_clock> startP,endP;
    startP = std::chrono::system_clock::now();
	test_hash(ab_arr, res, r);
    endP = std::chrono::system_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(endP-startP).count();
    std::cout << "duration hash TOTAL (seconds):  " << duration1/1e09 << std::endl;
    
    std::chrono::time_point<std::chrono::system_clock> startP2,endP2;
    startP2 = std::chrono::system_clock::now();
    test_array(ab_arr, res, r);
    endP2 = std::chrono::system_clock::now();
    duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(endP2-startP2).count();
    std::cout << "duration2 array TOTAL (seconds):  " << duration2/1e09 << std::endl;

    std::cout << std::endl;

    std::cout << std::numeric_limits<double>::digits10 << std::endl;
    std::cout << std::numeric_limits<float>::digits10 << std::endl;


	return 0;

}