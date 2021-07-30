#include "parlay/primitives.h"
#include <unordered_set>

template <class T>
parlay::sequence<T> dedup(parlay::sequence<T> const &_A) {
  auto A = _A;
  std::sort(A.begin(), A.end());
  A.erase(std::unique(A.begin(), A.end()), A.end());
  return A;
}
