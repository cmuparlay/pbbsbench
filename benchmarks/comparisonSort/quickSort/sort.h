#include "parlay/internal/quicksort.h"

constexpr bool INPLACE = true;

template <class T, class BinPred>
void compSort(parlay::sequence<T>  &A, const BinPred& f) {
  parlay::internal::p_quicksort_inplace(make_slice(A), f);
}
