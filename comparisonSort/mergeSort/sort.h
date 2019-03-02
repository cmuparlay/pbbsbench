#include "pbbslib/merge_sort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  pbbs::merge_sort_inplace(pbbs::range<E*>(A, A+n), f);
}
