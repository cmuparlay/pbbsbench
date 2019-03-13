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
#include <algorithm>
#include "parallel.h"
#include "geometry.h"
#include "sequence.h"
#include "hull.h"
using namespace std;

struct ipoint { point p; indexT idx;};

auto pntless = [&] (ipoint a, ipoint b) {
  return (a.p.x < b.p.x) || ((a.p.x == b.p.x) && (a.p.y < b.p.y));};

coord itriArea(ipoint a, ipoint b, ipoint c) {
  return triArea(a.p, b.p, c.p);}

pbbs::sequence<indexT> quickHull(pbbs::sequence<ipoint> P, ipoint l, ipoint r) {
  size_t n = P.size();
  //if (n <= 10000) return pbbs::sequence<indexT>();
  if (n <= 1) return pbbs::map<indexT>(P, [&] (ipoint p) {return p.idx;});
  else {
    size_t idx = max_element(pbbs::delayed_seq<coord>(n, [&] (size_t i) {
	  return itriArea(l ,r , P[i]);}), std::less<coord>());
    ipoint maxP = P[idx];

    pbbs::sequence<ipoint> left = pbbs::filter(P, [&] (ipoint p) {
	return itriArea(l, maxP, p) > 0;});
    
    pbbs::sequence<ipoint> right = pbbs::filter(P, [&] (ipoint p) {
      return itriArea(maxP, r, p) > 0;});
    P.clear();
    
    pbbs::sequence<indexT> leftR, rightR;
    par_do_if(n > 400,
	      [&] () {leftR = quickHull(std::move(left), l, maxP);},
	      [&] () {rightR = quickHull(std::move(right), maxP, r);},
	      true);
    
    pbbs::sequence<indexT> result(leftR.size() + rightR.size() + 1);
    pbbs::copy(leftR, result.slice(0, leftR.size()));
    result[leftR.size()] = maxP.idx;
    pbbs::copy(rightR, result.slice(leftR.size() + 1, result.size()));;
    return result;
  }
}

pbbs::sequence<indexT> hull(pbbs::sequence<point> const &Points) {
  timer t("hull", false);
  size_t n = Points.size();
  pbbs::sequence<ipoint> P(n, [&] (size_t i) {
      ipoint r = {Points[i], (indexT) i};
      return r;
    });

  t.next("convert");
  size_t mini, maxi;
  std::tie(mini, maxi) = pbbs::minmax_element(P, pntless);
  ipoint minp = P[mini];
  ipoint maxp = P[maxi];
  t.next("minmax");
  
  pbbs::sequence<bool> topFlag(n);
  pbbs::sequence<bool> bottomFlag(n);
  parallel_for (0, n, [&] (size_t i) {
    coord a = itriArea(minp, maxp, P[i]);
    topFlag[i] = a > 0;
    bottomFlag[i] = a < 0;
    });
  t.next("flags");

  pbbs::sequence<ipoint> top = pbbs::pack(P, topFlag);
  pbbs::sequence<ipoint> bottom = pbbs::pack(P, bottomFlag);
  P.clear();
  t.next("pack");

  pbbs::sequence<indexT> topR, bottomR;
  par_do([&] () {topR = quickHull(std::move(top), minp, maxp);},
	 [&] () {bottomR = quickHull(std::move(bottom), maxp, minp);}
	 );
  t.next("recurse");

  pbbs::sequence<indexT> result(topR.size() + bottomR.size() + 2);
  result[0] = minp.idx;
  pbbs::copy(topR, result.slice(1, 1 + topR.size()));
  result[1 + topR.size()] = maxp.idx;
  pbbs::copy(bottomR, result.slice(topR.size() + 2, result.size()));;
  t.next("append");
  return result;
}
