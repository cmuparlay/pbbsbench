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

#include "serialHull.h"

// The quickhull algorithm
// Points are all the points
// Idxs are the indices of the points (in Points) above the line defined by l-r.
// mid gives the index of the point furthest from the line defined by l-r
// The algorithm identifies the points above the lines l-mid and mid-r
//   and recurses on each
pbbs::sequence<indexT> quickHull(pbbs::sequence<point> const & Points,
				 pbbs::sequence<indexT> Idxs,
				 point l, indexT mid, point r) {
  size_t n = Idxs.size();
  if (n <= 1) return Idxs; 
  else {
    using cipair = std::pair<coord,indexT>;
    using cipairs = std::pair<cipair,cipair>;
    auto pairMax = [&] (cipairs a, cipairs b) {
      return cipairs((a.first.first > b.first.first) ? a.first : b.first,
		     (a.second.first > b.second.first) ? a.second : b.second);};

    // calculate furthest (positive) points from the lines l-mid and mid-r
    // at the same time set flags for those which are above each line
    pbbs::sequence<bool> leftFlag(n);
    pbbs::sequence<bool> rightFlag(n);
    point midP = Points[mid];
    auto P = pbbs::delayed_seq<cipairs>(n, [&] (size_t i) {
	indexT j = Idxs[i];
	coord lefta = triArea(l, midP, Points[j]);
	coord righta = triArea(midP, r, Points[j]);
	leftFlag[i] = lefta > 0.0;
	rightFlag[i] = righta > 0.0;
	return cipairs(cipair(lefta,j),cipair(righta,j));
      });
    cipairs prs = pbbs::reduce(P, pbbs::make_monoid(pairMax,cipairs()));
    indexT maxleft = prs.first.second;
    indexT maxright = prs.second.second;

    // keep those above each line
    pbbs::sequence<indexT> left = pbbs::pack(Idxs, leftFlag);
    pbbs::sequence<indexT> right = pbbs::pack(Idxs, rightFlag);
    Idxs.clear(); // clear and use std::move to avoid O(n log n) memory usage

    // recurse in parallel
    pbbs::sequence<indexT> leftR, rightR;
    par_do_if(n > 400,
	      [&] () {leftR = quickHull(Points, std::move(left), l, maxleft, midP);},
	      [&] () {rightR = quickHull(Points, std::move(right), midP, maxright, r);});

    // append the results together with mid in the middle
    pbbs::sequence<indexT> result(leftR.size() + rightR.size() + 1);
    pbbs::copy(leftR, result.slice(0, leftR.size()));
    result[leftR.size()] = mid;
    pbbs::copy(rightR, result.slice(leftR.size() + 1, result.size()));;
    return result;
  }
}

// The top-level call has to find the maximum and minimum x coordinates
//   and use them for the initial lines minp-maxp (for the upper hull)
//   and maxp-minp (for the lower hull).
pbbs::sequence<indexT> hull(pbbs::sequence<point> const &Points) {
  timer t("hull", true);
  size_t n = Points.size();
  auto pntless = [&] (point a, point b) {
    return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y));};

  // min and max points by x coordinate
  size_t minidx, maxidx, maxUpper, maxLower;
  std::tie(minidx, maxidx) = pbbs::minmax_element(Points, pntless);
  point minp = Points[minidx];
  point maxp = Points[maxidx];
  t.next("minmax");

  // identify those above and below the line minp-maxp
  // and calculate furtherst in each direction
  pbbs::sequence<bool> upperFlag(n);
  pbbs::sequence<bool> lowerFlag(n);
  auto P = pbbs::delayed_seq<coord>(n, [&] (size_t i) {
    coord a = triArea(minp, maxp, Points[i]);
    upperFlag[i] = a > 0;
    lowerFlag[i] = a < 0;
    return a;
    });
  std::tie(maxLower, maxUpper) = pbbs::minmax_element(P, std::less<coord>());
  
  t.next("flags");

  // pack the indices of those above and below
  pbbs::sequence<indexT> upper = pbbs::pack_index<indexT>(upperFlag);
  pbbs::sequence<indexT> lower = pbbs::pack_index<indexT>(lowerFlag);
  t.next("pack");

  // make parallel calls for upper and lower hulls
  pbbs::sequence<indexT> upperR, lowerR;
  par_do([&] () {upperR = quickHull(Points, std::move(upper), minp, maxUpper, maxp);},
	 [&] () {lowerR = quickHull(Points, std::move(lower), maxp, maxLower, minp);}
	 );
  t.next("recurse");

  pbbs::sequence<indexT> result(upperR.size() + lowerR.size() + 2);
  result[0] = minidx;
  pbbs::copy(upperR, result.slice(1, 1 + upperR.size()));
  result[1 + upperR.size()] = maxidx;
  pbbs::copy(lowerR, result.slice(upperR.size() + 2, result.size()));;
  t.next("append");
  return result;
}
