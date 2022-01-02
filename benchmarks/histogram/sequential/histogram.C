#include "parlay/primitives.h"
#include "parlay/io.h"

parlay::sequence<uint> histogram(parlay::sequence<uint> const &In, uint buckets) {
  parlay::sequence<uint> result(buckets+1);
  for (const auto& x : In) result[x]++;
  return result;
}
