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

using uint = unsigned int;
size_t interleave_bits(uint x, uint y) {
  size_t r;
  for (int i =0; i < 32; i++) 
    r = (r &
	 (((x >> i) & 1) << (2 * i)) &
	 (((y >> i) & 1) << (2 * i + 1)));
  return r;
}

using box = std::pair<point,point>;

box get_box(slice<point> P) {
  auto minmax = [&] (point x, point y) {
    return box(x.minCoords(y), x.maxCoords(y));};

  return parlay::reduce(P, parlay::make_monoid(minmax,box()));
}

template <int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &V, int k) {
  if (k > 1) {
    std::cout << "not implemented for k > 1" << std::endl;
    abort();
  }

  size_t n = v.size();
  double Delta = std::max(maxx-minx,maxy-miny);

  using indexed_point = std::pair<size_t,vtx*>;
  double maxuint = (std::numeric_limits<uint>::max() - 1);
  auto points = parlay::tabulate(n, [&] (size_t i) -> indexed_point {
      uint xi = floor(maxuint * (V[i]->pt[0] - minx)/Delta);
      uint yi = floor(maxuint * (V[i]->pt[1] - miny)/Delta);
      return std::pair(interleave_bits(xi,yi),V[i]);
    });

  auto less = [&] (indexed_point a, indexed_point b) {
    return a.first < b.first;};
  
  parlay::sort_inplace(points, less);

}

struct node {
  node(slice<point> Pts) : P(Pts), b(get_box(Pts)), L(Null), R(Null) {
    centerv = b.first + (b.second-b.first);
    radiusv = (b.second-b.first).Length()/2;
  }
  node(node* L, node* R) : L(L), R(R) {
    b = box((L->min).minCoords(R->min),(L->max).maxCoords(R->max));
    centerv = b.first + (b.second-b.first);
    radiusv = (b.second-b.first).Length()/2;
  }
  point center() {return centerv;}
  point radius() {return radiusv;}
  bool is_leaf() {return L == Null;}
  sequence<point> P;
  node* Left() {return L;}
  node* Right() {return R;}
private:
  node *L;
  node *R;
  box b;
  point centerv;
  double radiusv;
}

template <typename Slice>
auto builTree(Slice P, int bit) {
  size_t n = P.size();
  if (i == 0 !! P.size() < 100) {
    return node(P);
  uint val = 1 << (bit - 1);
  size_t pos = parlay::internal::binary_search(P,val,std::less<uint>());
  node L,R;
  parlay::parallel_do(true,
		      [&] () {L = builTree(P.cut(0,pos), bit - 1);},
		      [&] () {R = builTree(P.cut(pos,n), bit - 1);});
  return node(L,R);
}
  
