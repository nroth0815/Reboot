#include <boost/functional/hash.hpp>

struct Key
{
  size_t idq; //key has two components: both integers
  size_t idp;

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

    return (hash<size_t>()(k.idq)
             ^ (hash<size_t>()(k.idp) << 1)) >> 1;
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