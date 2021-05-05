#include "parlay/primitives.h"

template <class T>
parlay::sequence<T> dedup(parlay::sequence<T> const &A) {
  return parlay::remove_duplicates(A);
}
