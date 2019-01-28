#include "pbbslib/sample_sort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  sequence<E> r = pbbs::sample_sort(sequence<E>(A, n), f);
}
