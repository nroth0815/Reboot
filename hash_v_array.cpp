#include <string>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <boost/functional/hash.hpp> //better hash libraries (for later)
#include "hash_temp.hpp" //contains Key and Keyhasher*

#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

struct ab{

    MyFloat alpha;
    MyFloat beta;

} values;

void calc_ab(size_t res, ab *array){

    for(size_t i=0; i < res*res*res; i++){
        array[i].alpha=i;
        array[i].beta=-i;
    }

}

void test_hash(ab *array, size_t res){

    typedef std::unordered_map< Key, ab , KeyHasher_mod> hashmap;
    hashmap numbers;

    for(size_t i=0; i<res*res*res; i++){
        numbers[{i,i}]=array[i];
    }

    //numbers[{0,0}]=1;
    //numbers[{1,0}]=2;
    //numbers[{1243,230}]=-2;

    //std::cout << numbers[{0,0}]<< std::endl;

    hashmap::hasher hashfunc= numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
         std::cout << i->first.idq << ", " << i->first.idp << " -> "<< i->second.alpha << " " << second.beta << "; (hash = " << hashfunc( i->first ) << ")"  << std::endl;
    }


}

void test_array(ab *array, size_t res){


}

int main(int argc, char *argv[]){

    if(argc!=2){std::cerr << "Only one argument expected!" << std::endl; exit(1);}

    size_t res=atoi(argv[1]);
    
    ab *ab_arr=(ab*)calloc(res*res*res,sizeof(ab));

    //int ret=-1;
    calc_ab(res, ab_arr);

    uint64_t duration1, duration2;

    std::chrono::time_point<std::chrono::system_clock> startP,endP;
    startP = std::chrono::system_clock::now();
	test_hash();
    endP = std::chrono::system_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(endP-startP).count();
    std::cout << "duration1 (seconds):  " << duration1/1e09 << std::endl;
    
    std::chrono::time_point<std::chrono::system_clock> startP2,endP2;
    startP2 = std::chrono::system_clock::now();
    test_array();
    endP2 = std::chrono::system_clock::now();
    duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(endP2-startP2).count();
    std::cout << "duration2 (seconds):  " << duration2/1e09 << std::endl;

	return 0;

}