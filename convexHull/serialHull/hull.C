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
#include "geometry.h"
#include "hull.h"
#include "sequence.h"
using namespace std;

// bulk of code in subfile, so can be shared by parallel version
#include "serialHull.h"

pbbs::sequence<indexT> hull(pbbs::sequence<point> const &S) {
  point* P = S.begin();
  size_t n = S.size();
  indexT* I = pbbs::new_array<indexT>(n);
  for (size_t i=0; i < n; i++) I[i] = i;

  size_t l = 0;
  size_t r = 0;
  for (size_t i=1; i < n; i++) {
    if (P[i].x > P[r].x) r = i;
    if (P[i].x < P[l].x || ((P[i].x == P[l].x) && P[i].y < P[l].y)) 
      l = i;
  }

  pair<size_t,size_t> nn = split(I, n, aboveLine(P, l, r), aboveLine(P, r, l));
  size_t n1 = nn.first;
  size_t n2 = nn.second;

  size_t m1 = serialQuickHull(I, P, n1, l, r);
  size_t m2 = serialQuickHull(I+n-n2, P, n2, r, l);

  pbbs::sequence<indexT> RI(m1+2+m2);
  for (size_t i=m1; i > 0; i--) RI[i] = I[i-1];
  for (size_t i=0; i < m2; i++) RI[i+m1+2] = I[i+n-n2];
  RI[0] = l;
  RI[m1+1] = r;
  pbbs::free_array(I);
  return RI;
}

