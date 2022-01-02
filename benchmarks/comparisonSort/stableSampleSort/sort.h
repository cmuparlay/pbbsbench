#include "parlay/internal/sample_sort.h"

constexpr bool INPLACE = false;

template <class T, class BinPred>
parlay::sequence<T> compSort(parlay::sequence<T> const &A, const BinPred& f) {
  return parlay::internal::sample_sort(parlay::make_slice(A), f, true); // true makes it stable
}
