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

template <class ET, class F>
pair<size_t,size_t> split(ET* A, size_t n, F lf, F rf) {
  size_t ll = 0, lm = 0;
  size_t rm = n-1, rr = n-1;
  while (1) {
    while ((lm <= rm) && !(rf(A[lm]) > 0)) {
      if (lf(A[lm]) > 0) A[ll++] = A[lm];
      lm++;
    }
    while ((rm >= lm) && !(lf(A[rm]) > 0)) {
      if (rf(A[rm]) > 0) A[rr--] = A[rm];
      rm--;
    }
    if (lm >= rm) break; 
    ET tmp = A[lm++];
    A[ll++] = A[rm--];
    A[rr--] = tmp;
  }
  size_t n1 = ll;
  size_t n2 = n-rr-1;
  return pair<size_t,size_t>(n1,n2);
}

struct aboveLine {
  size_t l, r;
  point* P;
  aboveLine(point* _P, size_t _l, size_t _r) : P(_P), l(_l), r(_r) {}
  bool operator() (size_t i) {return triArea(P[l], P[r], P[i]) > 0.0;}
};

size_t serialQuickHull(indexT* I, point* P, size_t n, size_t l, size_t r) {
  if (n < 2) return n;
  size_t maxP = I[0];
  coord maxArea = triArea(P[l],P[r],P[maxP]);
  for (size_t i=1; i < n; i++) {
    size_t j = I[i];
    coord a = triArea(P[l],P[r],P[j]);
    if (a > maxArea) {
      maxArea = a;
      maxP = j;
    }
  }

  pair<size_t,size_t> nn = split(I, n, aboveLine(P,l,maxP), aboveLine(P,maxP,r));
  size_t n1 = nn.first;
  size_t n2 = nn.second;

  size_t m1, m2;
  m1 = serialQuickHull(I,      P, n1, l,   maxP);
  m2 = serialQuickHull(I+n-n2, P, n2, maxP,r);
  for (size_t i=0; i < m2; i++) I[i+m1+1] = I[i+n-n2];
  I[m1] = maxP;
  return m1+1+m2;
}

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

