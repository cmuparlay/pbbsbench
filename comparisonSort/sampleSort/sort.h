#include "parlay/internal/sample_sort.h"

<<<<<<< HEAD
constexpr bool INPLACE = false;

template <class T, class BinPred>
parlay::sequence<T> compSort(parlay::sequence<T> const &A, const BinPred& f) {
  return parlay::internal::sample_sort(parlay::make_slice(A), f);
=======
template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  parlay::internal::sample_sort(parlay::make_slice(A, A+ n), f);
>>>>>>> Magdalen
}
