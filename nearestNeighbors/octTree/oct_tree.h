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
#include "common/get_time.h"

// vtx must support v->pt
// and v->pt must support pt.dimension(), pt[i],
//    (pt1 - pt2).Length(), pt1 + (pt2 - pt3)
//    pt1.minCoords(pt2), pt1.maxCoords(pt2),
template <typename vtx>
struct oct_tree {

  using point = typename vtx::pointT;
  using uint = unsigned int;
  using box = std::pair<point,point>;
  using indexed_point = std::pair<size_t,vtx*>;
  using slice_t = decltype(make_slice(parlay::sequence<indexed_point>()));
  using slice_v = decltype(make_slice(parlay::sequence<vtx*>()));

  point center() {return centerv;}
  point radius() {return radiusv;}
  box Box() {return b;}
  size_t size() {return n;}
  bool is_leaf() {return L == nullptr;}
  parlay::sequence<vtx*> P;
  oct_tree* Left() {return L;}
  oct_tree* Right() {return R;}
  oct_tree* Parent() {return parent;}

private:
  constexpr static int key_bits = 64;
  double cut_plane;
  int cut_dim;
  size_t n;
  oct_tree *parent;
  oct_tree *L;
  oct_tree *R;
  box b;
  point centerv;
  double radiusv;
  static parlay::type_allocator<oct_tree> node_alloc;
    
  // takes a point, rounds each coordinate to an integer, and interleaves
  // the bits into "key_bits" total bits.
  // minp is the minimmum x,y,z coordinate for all points
  // delta is the largest range across the three dimensions
  static size_t interleave_bits(point p, point minp, double delta) {
    int dim = p.dimension();
    int bits = key_bits/dim;
    uint maxval = (((size_t) 1) << bits) - 1;
    uint ip[dim];
    for (int i = 0; i < dim; i++) 
      ip[i] = floor(maxval * (p[i] - minp[i])/delta);

    size_t r = 0;
    int loc = 0;
    for (int i =0; i < bits; i++)
      for (int d = 0; d < dim; d++) 
	r = r | (((ip[d] >> i) & (size_t) 1) << (loc++));
    //std::cout << "r = " << r << ", " << ip[0] << ", " << ip[1] <<std::endl;
    return r;
  }

  // generates a box consisting of a lower left corner,
  // and an upper right corner.
  static box get_box(parlay::sequence<vtx*> &V) {
    size_t n = V.size();
    auto minmax = [&] (box x, box y) {
      return box(x.first.minCoords(y.first),
		 x.second.maxCoords(y.second));};
    auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
	return box(V[i]->pt, V[i]->pt);});
    box identity = pts[0];
    return parlay::reduce(pts, parlay::make_monoid(minmax,identity));
  }

  static parlay::sequence<indexed_point> tag_points(parlay::sequence<vtx*> &V) {
    timer t("tag",false);
    size_t n = V.size();
    int dims = (V[0]->pt).dimension();
    
    box b = get_box(V);
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]);
    
    auto points = parlay::tabulate(n, [&] (size_t i) -> indexed_point {
	return std::pair(interleave_bits(V[i]->pt, b.first, Delta), V[i]);
      });
    t.next("tabulate");
    
    auto less = [] (indexed_point a, indexed_point b) {
      return a.first < b.first;};
    
    auto x = parlay::sort(points, less);
    t.next("sort");
    return x;
  }

  void set_center_radius() {			   
    centerv = b.first + (b.second-b.first)/2;
    //radiusv = (b.second-b.first).Length()/2;
  }

public:

  template <typename Seq>
  static oct_tree<vtx>* build(Seq &P) {
    timer t("oct_tree",false);
    int dims = (P[0]->pt).dimension();
    auto pts = tag_points(P);
    t.next("tag");
    //for (int i=0; i < P.size(); i++) 
    //  std::cout << pts[i].first << ", " << pts[i].second->pt[0] << ", ";
    //std::cout << std::endl;
    //std::cout << "here" << std::endl;
    oct_tree* r = node_alloc.alloc();
    new (r) oct_tree<vtx>(make_slice(pts), nullptr, dims*(key_bits/dims));
    t.next("build");
    return r;
  }

  oct_tree(slice_t Pts, oct_tree* parentp, int bit) {
    n = Pts.size();
    parent = parentp;
    int cutoff = 20;
    //std::cout << "size = " << n << std::endl;
    if (bit == 0 || Pts.size() < cutoff) {
      P = parlay::tabulate(Pts.size(), [&] (size_t i) -> vtx* {
	  return Pts[i].second;});
      L = R = nullptr;
      b = get_box(P);
    } else {
      size_t val = ((size_t) 1) << (bit - 1);
      size_t mask = (bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << bit);
      auto less = [&] (indexed_point x) {
	return (x.first & mask) < val;
      };
      size_t pos = parlay::internal::binary_search(Pts, less);
      // for (int i=0; i < n; i++) 
      //   std::cout << Pts[i].first << ", ";
      // std::cout << std::endl;
      // std::cout << "size = " << n << ", mask = " << mask << ", bit = "
      //           << bit << ", val = " << val << ", pos = " << pos << std::endl;
      if (pos == 0 || pos == n) {
	*this = oct_tree(Pts, parentp, bit - 1);
	return;
      }
      parlay::par_do(
		     ([&] () {
		       oct_tree* l = node_alloc.alloc();
		       L = new (l) oct_tree(Pts.cut(0,pos), this, bit - 1);}),
		     ([&] () {
		       oct_tree* r = node_alloc.alloc();
		       R = new (r) oct_tree(Pts.cut(pos,n), this, bit - 1);})
		     );
      b = box(L->b.first.minCoords(R->b.first),L->b.second.maxCoords(R->b.second));
    }
    set_center_radius();
    //std::cout << center()[0] << ", " << center()[1] << std::endl;
  }

  static void flatten_rec(oct_tree<vtx> *T, slice_v R) {
    if (T->is_leaf())
      for (int i=0; i < T->size(); i++)
	R[i] = T->P[i];
    else {
      size_t n_left = T->L->size();
      size_t n = T->size();
      parlay::par_do([&] () {flatten_rec(T->L, R.cut(0, n_left));},
		     [&] () {flatten_rec(T->R, R.cut(n_left, n));});
    }
  }

  template <typename F>
  void map(F f) {
    if (is_leaf())
      for (int i=0; i < size(); i++) f(P[i],this);
    else {
      parlay::par_do([&] () {L->map(f);},
		     [&] () {R->map(f);});
    }
  }

  size_t depth() {
    if (is_leaf()) return 0;
    else {
      size_t l, r;
      parlay::par_do([&] () {l = L->depth();},
		     [&] () {r = R->depth();});
      return 1 + std::max(l,r);
    }
  }
  
  parlay::sequence<vtx*> flatten() {
    parlay::sequence<vtx*> r(n);
    flatten_rec(this, parlay::make_slice(r));
    return r;
  }
};
