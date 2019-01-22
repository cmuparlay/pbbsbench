#pragma once
#include <parallel/algorithm>

template <class E, class BinPred>
void compSort(E* A, int n, BinPred f) {
  __gnu_parallel::sort(A,A+n,f);}
