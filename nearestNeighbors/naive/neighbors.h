// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"

int algorithm_version = 0; //necessary to make everything compile since octTree/neighbors.h requires it

// naive n^2 solution, currently just for k = 1
template <int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  if (k > 1) {
    std::cout << "not implemented for k > 1" << std::endl;
    abort();
  }
  size_t n = v.size();
  parlay::parallel_for (0, n, [&] (size_t i) {
      int k = (i + 1) % n;
      double d = (v[k]->pt - v[i]->pt).Length();
      for (int j = 0; j < n; j++) {
	if (j != i) 
	  if ((v[j]->pt - v[i]->pt).Length() < d) {
	    k = j;
	    d = (v[k]->pt - v[i]->pt).Length();
	  }
      }
      v[i]->ngh[0] = v[k];
    });
}
