#include <unordered_map>
#include <string>
#include <iostream>
#include <boost/functional/hash.hpp> //better hash libraries (for later)

//using namespace std;

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

	//typedef tr1::unordered_map< pair<int, int>, int, pair_hash> hashmap;
    typedef std::unordered_map< std::string, int> hashmap;
    hashmap numbers;

    //numbers[0]=1;
    //numbers[2]=-1;

    numbers["one"] = 1;
    numbers["two"] = 2;
    numbers["three"] = 3;

    //numbers[{0,0}]=1;
    //numbers[{1,0}]=2;

    //tr1::hash< int > hashfunc = numbers.hash_function();
    std::hash< std::string > hashfunc = numbers.hash_function();
    //tr1::hash< pair<int, int> > hashfunc = numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
        //cout << i << endl;
        std::cout << i->first << " -> " << i->second << " (hash = " << hashfunc( i->first ) << ")" << std::endl;
        //cout << i->first << " -> " << i->second << endl;//" (hash = " << pair_hash( i->first ) << ")" << endl;
    }

}

struct Key
{
  int first; //key has two components: both integers
  int second;

  bool operator==(const Key &other) const //self-defined comparison operator for our structure
  { return (first == other.first
            && second == other.second);
  }
};

namespace std { //adds the following to the standard namespace?

  template <>
  struct hash<Key>
  {
    std::size_t operator()(const Key& k) const
    {
      using std::size_t;
      using std::hash;
      //using std::int;

      // Compute individual hash values for first and
      // second and combine them using XOR
      // and bit shifting: (to avoid collisions)

      return (hash<int>()(k.first)
               ^ (hash<int>()(k.second) << 1)) >> 1;
    }
  };

}

struct KeyHasher
{
  std::size_t operator()(const Key& k) const
  {
    using std::size_t;
    using std::hash;
    //using std::string;

    return (hash<int>()(k.first)
             ^ (hash<int>()(k.second) << 1)) >> 1;
  }
};

void test_2d(){

    //typedef tr1::unordered_map< pair<int, int>, int, pair_hash> hashmap;
    typedef std::unordered_map< Key, int , KeyHasher> hashmap;
    hashmap numbers;

    //numbers[0]=1;
    //numbers[2]=-1;

    //numbers["one"] = 1;
    //numbers["two"] = 2;
    //numbers["three"] = 3;

    numbers[{0,0}]=1;
    numbers[{1,0}]=2;

    std::cout << numbers[{0,0}]<< std::endl;

    //tr1::hash< int > hashfunc = numbers.hash_function();
    // tr1::hash< string > hashfunc = numbers.hash_function();
    //tr1::hash< Key > hashfunc = numbers.hash_function();
    for( hashmap::const_iterator i = numbers.begin(), e = numbers.end() ; i != e ; ++i ) {
         std::cout << i->first.first << " " << i->first.second << std::endl;//" (hash = " << hashfunc( i->first ) << ")"  << endl;
    //     cout << i->first << " -> " << i->second << " (hash = " << hashfunc( i->first ) << ")" << endl;
    //     //cout << i->first << " -> " << i->second << endl;//" (hash = " << pair_hash( i->first ) << ")" << endl;
     }


     std::cout << "now to find items: "<< std::endl;
    //find items:
    //auto: deduces variable type from its initialiser, in this case: an iterator (new in c++11 standard)
    auto iter= numbers.find({0,0}); //this will throw an exception if key doesn't exist
    std::cout << iter->first.first << " " << iter->first.second << " " <<  iter->second << std::endl;

    //auto outp = numbers.at({1,0}); //returns the value directly; this will throw an exception if key doesn't exist
    //std::cout << outp << std::endl;



}

int main(){


	test_map();

    test_2d();

	return 0;


}