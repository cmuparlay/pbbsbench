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

#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"      
#include "parlay/primitives.h"
#include "parlay/alloc.h"
#include "common/geometry.h"
#include "common/get_time.h"
#include "flock/flock.h"
#include "verlib/verlib.h"

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

  constexpr static int node_cutoff = 32;
  constexpr static int key_bits = 64;

  static size_t interleave_bits(point p, point min_point, double delta) {
    int dim = p.dimension();
    int bits = key_bits/dim;
    uint maxval = (((size_t) 1) << bits) - 1; //maybe should just be size_t instead of uint
    uint ip[dim];
    for (int i = 0; i < dim; i++) 
      ip[i] = floor(maxval * (p[i] - min_point[i])/delta); //could be something other than floor? nearest to?
    size_t r = 0;
    int loc = 0;
    for (int i =0; i < bits; i++)
      for (int d = 0; d < dim; d++) 
  r = r | (((ip[d] >> i) & (size_t) 1) << (loc++));
    return r;
  }

  // generates a box consisting of a lower left corner,
  // and an upper right corner.
  static box get_box(parlay::sequence<vtx*> &V) { // parlay::sequence<vtx*> &V) {
    size_t n = V.size();
    auto minmax = [&] (box x, box y) {
      return box(x.first.minCoords(y.first),
     x.second.maxCoords(y.second));};
    // uses a delayed sequence to avoid making a copy
    auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
      return box(V[i]->pt, V[i]->pt);});
    box identity = pts[0];  
    box final = parlay::reduce(pts, parlay::make_monoid(minmax, identity));
    return (final);
  }

    // generates a box consisting of a lower left corner,
  // and an upper right corner.
  static box get_box(parlay::sequence<indexed_point> &V) { // parlay::sequence<vtx*> &V) {
    size_t n = V.size();
    // std::cout << n << std::endl;
    auto minmax = [&] (box x, box y) {
      return box(x.first.minCoords(y.first),
     x.second.maxCoords(y.second));};
    // uses a delayed sequence to avoid making a copy
    auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
      return box(V[i].second->pt, V[i].second->pt);});
    box identity = pts[0];  
    box final = parlay::reduce(pts, parlay::make_monoid(minmax, identity));
    // node::print_point(final.first);
    // node::print_point(final.second);
    return (final);
  }

  struct node : verlib::versioned { 

  public:
    int bit;
    flck::atomic_write_once<bool> removed;
    flck::lock lck;
    // using leaf_seq = parlay::sequence<vtx*>;
    point center() {return centerv;}
    box Box() {return b;}
    size_t size() {return n;} //NOT ALWAYS ACCURATE
    bool is_leaf() {return (L.load() == nullptr) && (R.load() == nullptr);}
    node* Left() {return L.load();}
    node* Right() {return R.load();}
    parlay::sequence<indexed_point>& Indexed_Pts(){return indexed_pts;}
    // leaf_seq& Vertices() {return P;}

    void set_size(size_t s){n=s;}

    void set_child(node* child, bool left){
      if(left) L = child;
      else R = child;
    }

    static void print_point(point p){
      int d = p.dimension();
      std::cout << "Point: ";
      for(int j=0; j<d; j++){
        std::cout << p[j] << ", ";
      }
      std::cout << "\n";
    }

    // construct a leaf node with a sequence of points directly in it
    node(parlay::sequence<indexed_point> &Pts, int currentBit) : removed(false), n(Pts.size()) { 
      // strips off the integer tag, no longer needed
      indexed_pts = std::move(Pts);
      L = R = nullptr;
      b = get_box(indexed_pts);
      set_center();
      bit = currentBit;
    }

    // construct a leaf node with a sequence of points directly in it
    node(parlay::sequence<indexed_point> &Pts, int currentBit, box &B) : removed(false), n(Pts.size()) { 
      // strips off the integer tag, no longer needed
      indexed_pts = std::move(Pts);
      L = R = nullptr;
      b = B;
      set_center();
      bit = currentBit;
    }

    // construct a leaf node with a sequence of points directly in it
    node(slice_t Pts, int currentBit, box &B) : removed(false), n(Pts.size()) { 
      // strips off the integer tag, no longer needed
      indexed_pts = parlay::sequence<indexed_point>(n);
      for (int i = 0; i < n; i++) {
        indexed_pts[i] = Pts[i];  
      }
      L = R = nullptr;
      b = B;
      set_center();
      bit = currentBit;
    }


    // construct a leaf node with a sequence of points directly in it
    node(slice_t Pts, int currentBit) : removed(false), n(Pts.size()) { 
      // strips off the integer tag, no longer needed
      indexed_pts = parlay::sequence<indexed_point>(n);
      for (int i = 0; i < n; i++) {
        if(indexed_pts.size() < n) {
          std::stringstream ss;
          ss << "n = " << n << std::endl << "P.size() = " << indexed_pts.size() << std::endl;
          std::cout << ss.str() << std::endl;
          exit(1);
        }
        indexed_pts[i] = Pts[i];  
      }
      L = R = nullptr;
      b = get_box(indexed_pts);
      set_center();
      bit = currentBit;
    }

    // construct an internal binary node
    node(node* L, node* R, int currentBit) : removed(false), L(L), R(R) { 
      b = box(L->b.first.minCoords(R->b.first),
	      L->b.second.maxCoords(R->b.second));
      n = L->size() + R->size();
      set_center();
      bit = currentBit;
    }

    node(node* L, node* R, int currentBit, box B) : removed(false), L(L), R(R) { 
      b = B;
      n = L->size() + R->size();
      set_center();
      bit = currentBit;
    }

    // recursively frees the tree
    static void delete_rec(node* T) {
      if (T != nullptr) {
        parlay::par_do_if(T->n > 1000,
         [&] () { delete_rec(T->L.load());},
         [&] () { delete_rec(T->R.load());});
        node::free_node(T);
      }
    }

    static node* new_leaf(slice_t Pts, int currentBit) {
      node* r = alloc_node();
      assert(Pts.begin() != nullptr);
      new (r) node(Pts, currentBit);
      return r;
    }

    static node* new_leaf(slice_t Pts, int currentBit, box B) {
      node* r = alloc_node();
      assert(Pts.begin() != nullptr);
      new (r) node(Pts, currentBit, B);
      return r;
    }

    static node* new_leaf(parlay::sequence<indexed_point> Pts, int currentBit) {
      node* r = alloc_node();
      assert(Pts.begin() != nullptr);
      new (r) node(Pts, currentBit);
      return r;
    }

    static node* new_leaf(parlay::sequence<indexed_point> Pts, int currentBit, box B) {
      node* r = alloc_node();
      assert(Pts.begin() != nullptr);
      new (r) node(Pts, currentBit, B);
      return r;
    }

    static node* new_node(node* L, node* R, int currentBit) {
      node* nd = alloc_node();
      new (nd) node(L, R, currentBit);
      return nd;
    }

    static node* new_node(node* L, node* R, int currentBit, box B) {
      node* nd = alloc_node();
      new (nd) node(L, R, currentBit, B);
      return nd;
    }
    
    // ~node() {
    //   // need to collect in parallel
    //   parlay::par_do_if(n > 1000,
		// 	[&] () { delete_tree(L.load());},
		// 	[&] () { delete_tree(R.load());});
    // }

    parlay::sequence<vtx*> flatten() {
      parlay::sequence<vtx*> r(n);
      flatten_rec(this, parlay::make_slice(r));
      return r;
    }

    // map a function f(p,node_ptr) over the points, passing
    // in a pointer to a vertex, and a pointer to the leaf node it is in.
    // f should return void

    //pass in a function to compute nearest neighbors
    template <typename F>
    void map(F f) { 
      if (is_leaf())
	for (int i=0; i < size(); i++) f(indexed_pts[i].second,this);
      else {
	parlay::par_do_if(n > 1000,
			  [&] () {L.load()->map(f);},
			  [&] () {R.load()->map(f);});
      }
    }


    size_t depth() {
      if (is_leaf()) return 0;
      else {
	size_t l, r;
	parlay::par_do_if(n > 1000,
			  [&] () {l = L.load()->depth();},
			  [&] () {r = R.load()->depth();});
	return 1 + std::max(l,r);
      }
    }

    // disable copy and move constructors/assignment since
    // they are dangerous with with free.
    node(const node&) = delete;
    node(node&&) = delete;
    node& operator=(node const&) = delete;
    node& operator=(node&&) = delete;
    node() { }

    static node* alloc_node();
    static void free_node(node* T);
    static void retire_node(node* T);
    static void shuffle(size_t n);
  private:

    size_t n; //NOT ALWAYS ACCURATE
    verlib::versioned_ptr<node> L;
    verlib::versioned_ptr<node> R;
    box b;
    point centerv;
    parlay::sequence<indexed_point> indexed_pts;
    // leaf_seq P;

    void set_center() {			   
      centerv = b.first + (b.second-b.first)/2;
    }

    static void flatten_rec(node *T, slice_v R) {
      if (T->is_leaf())
	for (int i=0; i < T->size(); i++)
	  R[i] = T->P[i];
      else {
	size_t n_left = T->L.load()->size();
	size_t n = T->size();
	parlay::par_do_if(n > 1000,
	  [&] () {flatten_rec(T->L, R.cut(0, n_left));},
	  [&] () {flatten_rec(T->R, R.cut(n_left, n));});
      }
    }
  }; // end struct node

  // build a tree given a sequence of pointers to points
  template <typename Seq>
  static node* build(Seq &P) {
    timer t("oct_tree", false);
    int dims = (P[0]->pt).dimension();
    auto pts = tag_points(P);
    t.next("tag");
    node* r = build_recursive(make_slice(pts), dims*(key_bits/dims));
    t.next("build");
    return r;
  }

    // build a tree given a sequence of pointers to points
  template <typename Seq>
  static node* build(Seq &P, box b) {
    timer t("oct_tree", false);
    int dims = (P[0]->pt).dimension();
    auto pts = tag_points(P, b);
    t.next("tag");
    node* r = build_recursive(make_slice(pts), dims*(key_bits/dims));
    t.next("build");
    return r;
  }

private:
  
 

  // tags each point (actually a pointer to it), with an interger
  // consisting of the interleaved bits for the x,y,z coordinates.
  // Also sorts based the integer.
  static parlay::sequence<indexed_point> tag_points(parlay::sequence<vtx*> &V) {
    timer t("tag", false); //tag is an arbitrary string, turn to true for printing out
    size_t n = V.size();
    int dims = (V[0]->pt).dimension();

    // find box around points, and size along largest axis
    box b = get_box(V);
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]); 
    t.next("get box");
    
    auto points = parlay::delayed_seq<indexed_point>(n, [&] (size_t i) -> indexed_point {
	return std::pair(interleave_bits(V[i]->pt, b.first, Delta), V[i]);
      }); //make this not a delayed sequence, tabulate instead, so that we can use t.next()
    
    auto less = [] (indexed_point a, indexed_point b) {
      return a.first < b.first;};
    
    auto x = parlay::sort(points, less);
    t.next("tabulate and sort");
    return x;
  }

    static parlay::sequence<indexed_point> tag_points(parlay::sequence<vtx*> &V, box b) {
    timer t("tag", false); //tag is an arbitrary string, turn to true for printing out
    size_t n = V.size();
    int dims = (V[0]->pt).dimension();

    // find box around points, and size along largest axis
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]); 
    t.next("get box");
    
    auto points = parlay::delayed_seq<indexed_point>(n, [&] (size_t i) -> indexed_point {
  return std::pair(interleave_bits(V[i]->pt, b.first, Delta), V[i]);
      }); //make this not a delayed sequence, tabulate instead, so that we can use t.next()
    
    auto less = [] (indexed_point a, indexed_point b) {
      return a.first < b.first;};
    
    auto x = parlay::sort(points, less);
    t.next("tabulate and sort");
    return x;
  }

  // each point is a pair consisting of an interleave integer along with
  // the pointer to the point.   The bit specifies which bit of the integer
  // we are working on (starts at top, and goes down).
  static node* build_recursive(slice_t Pts, int bit) {
    size_t n = Pts.size();
    if (n == 0) abort();

    // if run out of bit, or small then generate a leaf
    if (bit == 0 || n < node_cutoff) {
      return node::new_leaf(Pts, bit); 
    } else {

      // this was extracted to lookup_bit but left as is here since the less function requires mask and val
      size_t val = ((size_t) 1) << (bit - 1);
      size_t mask = (bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << bit);
      auto less = [&] (indexed_point x) {
	return (x.first & mask) < val;
      };
      // and then we binary search for the cut point
      size_t pos = parlay::internal::binary_search(Pts, less);

      // if all points are on one side, then move onto the next bit
      if (pos == 0 || pos == n) 
	      return build_recursive(Pts, bit - 1);

      // otherwise recurse on the two parts, also moving to next bit
      else {
	node *L, *R; 
	parlay::par_do_if(n > 1000,
           [&] () {L = build_recursive(Pts.cut(0, pos), bit - 1);},
	   [&] () {R = build_recursive(Pts.cut(pos, n), bit - 1);});
	return node::new_node(L, R, bit); 
      }
    }
  }
}; //end octTree structure 

  // uses the parlay memory manager, could be replaced will alloc/free

template <typename vtx>
flck::memory_pool<typename oct_tree<vtx>::node> node_allocator;

template <typename vtx>
typename oct_tree<vtx>::node* oct_tree<vtx>::node::alloc_node() { return node_allocator<vtx>.new_obj();}

template <typename vtx>
void oct_tree<vtx>::node::free_node(node* T) { node_allocator<vtx>.destruct(T);}
  
template <typename vtx>
void oct_tree<vtx>::node::retire_node(node* T) { node_allocator<vtx>.retire(T);}

template <typename vtx>
void oct_tree<vtx>::node::shuffle(size_t n) { node_allocator<vtx>.shuffle(10000000);}


  
