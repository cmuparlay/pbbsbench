#include "pbbslib/quicksort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  pbbs::slice_t<E*> B(A,A+n);
  pbbs::p_quicksort_inplace(B, f);
}
