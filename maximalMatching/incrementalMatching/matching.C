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

#define NOTMAIN 1
#include <iostream>
#include "sequence.h"
#include "parallel.h"
#include "graph.h"
#include "speculative_for.h"
#include "get_time.h"
#include "matching.h"
using namespace std;

using reservation = pbbs::reservation<edgeId>;

struct matchStep {
  edges& E;
  pbbs::sequence<reservation> &R;
  pbbs::sequence<bool> &matched;

  matchStep(edges E,
	    pbbs::sequence<reservation> &R,
	    pbbs::sequence<bool> &matched)
    : E(E), R(R), matched(matched) {}

  bool reserve(edgeId i) {
    size_t u = E[i].u;
    size_t v = E[i].v;
    if (matched[u] || matched[v] || (u == v)) return 0;
    R[u].reserve(i);
    R[v].reserve(i);
    return 1;
  }

  bool commit(edgeId i) {
    size_t u = E[i].u;
    size_t v = E[i].v;
    if (R[v].check(i)) {
      R[v].reset();
      if (R[u].check(i)) {
	matched[u] = matched[v] = 1;
	return 1;
      } 
    } else if (R[u].check(i)) R[u].reset();
    return 0;
  }
};

pbbs::sequence<edgeId> maximalMatching(edges E) {
  size_t n = max(E.numCols,E.numRows);
  size_t m = E.nonZeros;

  pbbs::sequence<reservation> R(n);
  pbbs::sequence<bool> matched(n, false);
  matchStep mStep(E, R, matched);
  pbbs::speculative_for<edgeId>(mStep, 0, m, 50, 0);
  pbbs::sequence<edgeId> matchingIdx =
    pbbs::pack(pbbs::delayed_seq<edgeId>(n, [&] (size_t i) {return R[i].r;}),
	       pbbs::sequence<bool>(n, [&] (size_t i) {return R[i].reserved();}));
  cout << "number of matches = " << matchingIdx.size() << endl;
  return matchingIdx;
}  
