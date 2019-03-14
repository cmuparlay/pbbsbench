#include "hull.h"
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
