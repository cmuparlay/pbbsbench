#include "pbbslib/sample_sort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  pbbs::sample_sort(A, n, f);
}
