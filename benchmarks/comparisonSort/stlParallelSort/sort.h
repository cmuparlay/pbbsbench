#include <parallel/algorithm>
#include "parlay/sequence.h"

constexpr bool INPLACE = true;

template <class T, class BinPred>
void compSort(parlay::sequence<T> &A, const BinPred& f) {
  __gnu_parallel::sort(A.begin(), A.end(), f);}
