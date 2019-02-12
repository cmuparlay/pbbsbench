#include "pbbslib/merge_sort.h"
#include "pbbslib/seq.h"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  pbbs::sequence<E> In(A,n);
  pbbs::sequence<E> Out = pbbs::merge_sort(std::move(In), f);
}
