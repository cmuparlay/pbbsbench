#include "parlay/internal/merge_sort.h"

constexpr bool INPLACE = true;

template <class T, class BinPred>
void compSort(parlay::sequence<T> &A, const BinPred& f) {
  parlay::internal::merge_sort_inplace(parlay::make_slice(A), f);
}
