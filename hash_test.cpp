#include <unordered_map>
#include <string>
#include <iostream>
#include <boost/functional/hash.hpp> //better hash libraries (for later)
#include "hash_temp.hpp" //contains Key and Keyhasher*

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};

void test_map(){

	typedef std::unordered_map< std::string, int> hashmap;
    hashmap numbers;

    numbers["one"] = 1;
    numbers["two"] = 2;
    numbers["three"] = 3;

    std::hash< std::string > hashfunc = numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
        std::cout << i->first << " -> " << i->second << " (hash = " << hashfunc( i->first ) << ")" << std::endl;
    }
    std::cout<< std::endl;

}

// struct Key
// {
//   int idq; //key has two components: both integers
//   int idp;

//   bool operator==(const Key &other) const //self-defined comparison operator for our structure
//   { return (idq == other.idq
//             && idp == other.idp);
//   }
// };

// struct KeyHasher //the actual hash calculation
// {
//   std::size_t operator()(const Key& k) const
//   {
//     using std::size_t;
//     using std::hash;

//     return (hash<int>()(k.idq)
//              ^ (hash<int>()(k.idp) << 1)) >> 1;
//     // Mainly for demonstration purposes, i.e. works but is overly simple
//     // In the real world, use sth. like boost.hash_combine
//   }
// };

// struct KeyHasher_mod //the actual hash calculation
// {
//   std::size_t operator()(const Key& k) const
//   {
//     using std::size_t;
//     using std::hash;
//     std::size_t seed=0;

//     boost::hash_combine(seed, k.idq);
//     boost::hash_combine(seed, k.idp);

//     return seed;
  
//   }
// };

void test_2d(){

    typedef std::unordered_map< Key, int , KeyHasher> hashmap;
    hashmap numbers;

    numbers[{0,0}]=1;
    numbers[{1,0}]=2;
    numbers[{1243,230}]=-2;

    std::cout << numbers[{0,0}]<< std::endl;

    hashmap::hasher hashfunc= numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
         std::cout << i->first.idq << ", " << i->first.idp << " -> "<< i->second << "; (hash = " << hashfunc( i->first ) << ")"  << std::endl;
    }


    std::cout << "now to find items: "<< std::endl;
    //find items:
    //auto: deduces variable type from its initialiser, in this case: an iterator (new in c++11 standard)
    auto iter= numbers.find({0,0}); //this will throw an exception if key doesn't exist
    std::cout << iter->first.idq << " " << iter->first.idp << " " <<  iter->second << std::endl;

    auto outp = numbers.at({1,0}); //returns the value directly; this will throw an exception if key doesn't exist
    std::cout << outp << std::endl;
    std::cout<< std::endl;

}

void test_2d_mod(){

    typedef std::unordered_map< Key, int , KeyHasher_mod> hashmap;
    hashmap numbers;

    numbers[{0,0}]=1;
    numbers[{1,0}]=2;
    numbers[{1243,230}]=-2;

    std::cout << numbers[{0,0}]<< std::endl;

    hashmap::hasher hashfunc= numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
         std::cout << i->first.idq << ", " << i->first.idp << " -> "<< i->second << "; (hash = " << hashfunc( i->first ) << ")"  << std::endl;
    }


    std::cout << "now to find items: "<< std::endl;
    //find items:
    //auto: deduces variable type from its initialiser, in this case: an iterator (new in c++11 standard)
    auto iter= numbers.find({0,0}); //this will throw an exception if key doesn't exist
    std::cout << iter->first.idq << " " << iter->first.idp << " " <<  iter->second << std::endl;

    auto outp = numbers.at({1,0}); //returns the value directly; this will throw an exception if key doesn't exist
    std::cout << outp << std::endl;
    std::cout<< std::endl;

}

int main(){

	test_map();
    test_2d();
    test_2d_mod();

	return 0;

}