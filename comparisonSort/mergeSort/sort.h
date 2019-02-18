#include "pbbslib/merge_sort.h"
#include "pbbslib/seq.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  pbbs::merge_sort_inplace(pbbs::slice_t<E*>(A, A+n), f);
}
