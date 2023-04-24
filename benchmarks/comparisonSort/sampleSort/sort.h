#include "sample_sort.h"

constexpr bool INPLACE = true;

template <class T, class BinPred>
void compSort(parlay::sequence<T> &A, const BinPred& f) {
  sample_sort_inplace(A, f);
}

// The following also works but does not sort in place
// performance is about the same, but a little slower on strings due to extra copy.

// constexpr bool INPLACE = false;

// template <class T, class BinPred>
// auto compSort(parlay::sequence<T> &A, const BinPred& f) {
//   return sample_sort(A, f);
// }

