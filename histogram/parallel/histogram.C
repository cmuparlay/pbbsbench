#include "parlay/primitives.h"
#include "parlay/io.h"

parlay::sequence<uint> histogram(parlay::sequence<uint> const &In, uint buckets) {
  auto x = parlay::tabulate(In.size(), [&] (uint i) {return std::pair(In[i]%256,i);});
  return parlay::histogram(In, buckets);
}
