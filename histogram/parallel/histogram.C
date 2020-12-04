#include "parlay/primitives.h"

parlay::sequence<uint> histogram(parlay::sequence<uint> const &In, uint buckets) {
  return parlay::histogram(In, buckets);
}
