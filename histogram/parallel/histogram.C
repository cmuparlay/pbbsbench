#include "parlay/primitives.h"
#include "parlay/io.h"

parlay::sequence<uint> histogram(parlay::sequence<uint> const &In, uint buckets) {
  return parlay::histogram(In, buckets);
}
