#include <algorithm>

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  std::sort(A, A+ n, f);
}
