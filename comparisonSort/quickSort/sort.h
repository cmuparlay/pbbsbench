#include "parlay/internal/quicksort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  parlay::slice<E*,E*> B(A,A+n);
  parlay::internal::p_quicksort_inplace(B, f);
}
