#include "pbbslib/quicksort.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  E* B = pbbs::new_array_no_init<E>(n,1);
  pbbs::p_quicksort(sequence<E>(A,n), sequence<E>(B,n), f, 1);
  pbbs::free_array(B);
}
