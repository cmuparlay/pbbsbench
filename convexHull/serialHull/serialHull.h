#include "hull.h"
using namespace std;

template <class ET, class FL, class FR>
pair<long,long> split(ET* A, long n, FL lf, FR rf) {
  long ll = 0, lm = 0;
  long rm = n-1, rr = n-1; // need to be signed
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
  long n1 = ll;
  long n2 = n-rr-1;
  return pair<long,long>(n1,n2);
}

indexT serialQuickHull(indexT* I, const point* P, indexT n, indexT l, indexT r) {
  if (n < 2) return n;
  indexT maxP = I[0];
  coord maxArea = triArea(P[l],P[r],P[maxP]);
  for (indexT i=1; i < n; i++) {
    indexT j = I[i];
    coord a = triArea(P[l],P[r],P[j]);
    if (a > maxArea) {
      maxArea = a;
      maxP = j;
    }
  }

  auto aboveLeft = [&] (indexT i) {
    return triArea(P[l], P[maxP], P[i]) > 0.0;};
  auto aboveRight = [&] (indexT i) {
    return triArea(P[maxP], P[r], P[i]) > 0.0;};

  pair<indexT,indexT> nn = split(I, n, aboveLeft, aboveRight);
  indexT n1 = nn.first;
  indexT n2 = nn.second;

  indexT m1, m2;
  m1 = serialQuickHull(I,      P, n1, l,   maxP);
  m2 = serialQuickHull(I+n-n2, P, n2, maxP,r);
  for (indexT i=0; i < m2; i++) I[i+m1+1] = I[i+n-n2];
  I[m1] = maxP;
  return m1+1+m2;
}
