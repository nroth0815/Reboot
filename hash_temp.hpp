#include <boost/functional/hash.hpp>

struct Key
{
  int idq; //key has two components: both integers
  int idp;

  bool operator==(const Key &other) const //self-defined comparison operator for our structure
  { return (idq == other.idq
            && idp == other.idp);
  }
};

struct KeyHasher //the actual hash calculation
{
  std::size_t operator()(const Key& k) const
  {
    using std::size_t;
    using std::hash;

    return (hash<int>()(k.idq)
             ^ (hash<int>()(k.idp) << 1)) >> 1;
    // Mainly for demonstration purposes, i.e. works but is overly simple
    // In the real world, use sth. like boost.hash_combine
  }
};

struct KeyHasher_mod //the actual hash calculation
{
  std::size_t operator()(const Key& k) const
  {
    using std::size_t;
    using std::hash;
    std::size_t seed=0;

    boost::hash_combine(seed, k.idq);
    boost::hash_combine(seed, k.idp);

    return seed;
  
  }
};