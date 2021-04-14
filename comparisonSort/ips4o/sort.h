#include "ips4o.hpp"
#include "parlay/sequence.h"

template <class T, class BinPred>
void compSort(parlay::sequence<T> &A, const BinPred& f) {
  ips4o::parallel::sort(A.begin(), A.end(), f);
}
