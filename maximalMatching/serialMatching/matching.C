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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/graph.h"
#include "common/speculative_for.h"
#include "common/get_time.h"
#include "matching.h"
using namespace std;

parlay::sequence<edgeId> maximalMatching(edges const &E) {
  size_t n = max(E.numCols,E.numRows);
  size_t m = E.nonZeros;
  timer t("max matching", true);
  parlay::sequence<edgeId> matching(n);
  parlay::sequence<bool> matched(n, false);
  size_t offset = 0;
  for (size_t i=0; i<m; i++) {
    size_t u = E[i].u;
    size_t v = E[i].v;
    if (matched[u] || matched[v]) continue;
    matched[u] = true;
    matched[v] = true;
    matching[offset] = i;
    offset++;
  }
  matching.resize(offset);
  return matching;
}
