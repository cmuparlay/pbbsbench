#include "parlay/internal/merge_sort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  parlay::internal::merge_sort_inplace(parlay::make_slice(A, A+n), f);
}
