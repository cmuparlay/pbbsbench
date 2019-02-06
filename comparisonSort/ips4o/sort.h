#include "ips4o.hpp"

template <class E, class BinPred>
void compSort(E* A, unsigned int n, const BinPred& f) {
  ips4o::parallel::sort(A, A+n, f);
}
