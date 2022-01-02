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

#define report_stats false
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "oct_tree.h"

// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vtx, int max_k>
struct k_nearest_neighbors {
  using point = typename vtx::point_t;
  using fvect = typename point::vector;
  using o_tree = oct_tree<vtx>;
  using node = typename o_tree::node;
  using tree_ptr = typename o_tree::tree_ptr;

  tree_ptr tree;

  // generates the search structure
  k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    tree = o_tree::build(V);
  }

  // returns the vertices in the search structure, in an
  //  order that has spacial locality
  parlay::sequence<vtx*> vertices() {
    return tree->flatten();
  }

  struct kNN {
    vtx *vertex;  // the vertex for which we are trying to find a NN
    vtx *neighbors[max_k];  // the current k nearest neighbors (nearest last)
    double distances[max_k]; // distance to current k nearest neighbors
    int k;
    int dimensions;
    size_t leaf_cnt;
    size_t internal_cnt;
    kNN() {}

    // returns the ith smallest element (0 is smallest) up to k-1
    vtx* operator[] (const int i) { return neighbors[k-i-1]; }

    kNN(vtx *p, int kk) {
      if (kk > max_k) {
	std::cout << "k too large in kNN" << std::endl;
	abort();}
      k = kk;
      vertex = p;
      dimensions = p->pt.dimension();
      leaf_cnt = internal_cnt = 0;
      // initialize nearest neighbors to point to Null with
      // distance "infinity".
      for (int i=0; i<k; i++) {
	neighbors[i] = (vtx*) NULL; 
	distances[i] = numeric_limits<double>::max();
      }
    }

    // if p is closer than neighbors[0] then swap it in
    void update_nearest(vtx *other) { 
      auto dist = (vertex->pt - other->pt).Length();
      if (dist < distances[0]) {
	neighbors[0] = other;
	distances[0] = dist;
	for (int i = 1;
	     i < k && distances[i-1] < distances[i];
	     i++) {
	  swap(distances[i-1], distances[i]);
	  swap(neighbors[i-1], neighbors[i]); }
      }
    }

    bool within_epsilon_box(node* T, double epsilon) {
      auto box = T->Box();
      bool result = true;
      for (int i = 0; i < dimensions; i++) {
	result = (result &&
		  (box.first[i] - epsilon < vertex->pt[i]) &&
		  (box.second[i] + epsilon > vertex->pt[i]));
      }
      return result;
    }

    double distance(node* T) {
      return (T->center() - vertex->pt).Length();
    }
    
    // looks for nearest neighbors for this->vertex in Tree node T
    void k_nearest_rec(node* T) {
      if (report_stats) internal_cnt++;
      if (within_epsilon_box(T, distances[0])) {
	if (T->is_leaf()) {
	  if (report_stats) leaf_cnt++;
	  auto &Vtx = T->Vertices();
	  for (int i = 0; i < T->size(); i++)
	    if (Vtx[i] != vertex) update_nearest(Vtx[i]);
	} else if (distance(T->Left()) < distance(T->Right())) {
	  k_nearest_rec(T->Left());
	  k_nearest_rec(T->Right());
	} else {
	  k_nearest_rec(T->Right());
	  k_nearest_rec(T->Left());
	}
      }
    }

    // finds a point that is vaguely near
    void near_rec(node* T) {
      if (T->is_leaf()) {
	auto &Vtx = T->Vertices();
	for (int i = 0; i < T->size(); i++)
	  if (Vtx[i] != vertex) update_nearest(Vtx[i]);
      } else if (distance(T->Left()) < distance(T->Right())) {
	near_rec(T->Left());
      } else {
	near_rec(T->Right());
      }
    }

  };

  void k_nearest(vtx *p, int k) {
    kNN nn(p,k);
    nn.k_nearest_rec(tree.get());
    if (report_stats) p->counter = nn.internal_cnt;
    for (int i=0; i < k; i++)
      p->ngh[i] = nn[i];
  }
  
  vtx* nearest(vtx *p) {
    kNN nn(p,1);
    nn.k_nearest_rec(tree.get());
    if (report_stats) p->counter = nn.internal_cnt;
    return nn[0];
  }

  vtx* near(vtx *p) {
    kNN nn(p,1);
    nn.near_rec(tree.get());
    return nn[0];
  }

};

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int max_k, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN",report_stats);

  {
    using knn_tree = k_nearest_neighbors<vtx, max_k>;
    knn_tree T(v);
    t.next("build tree");

    if (report_stats) 
      std::cout << "depth = " << T.tree->depth() << std::endl;

    // this reorders the vertices for locality
    parlay::sequence<vtx*> vr = T.vertices();
    t.next("flatten tree");
  
    // find nearest k neighbors for each point
    parlay::parallel_for (0, v.size(), [&] (size_t i) {
					 T.k_nearest(vr[i], k);}, 1);

    t.next("try all");
    if (report_stats) {
      auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
      size_t i = parlay::max_element(s) - s.begin();
      size_t sum = parlay::reduce(s);
      std::cout << "max internal = " << s[i] 
		<< ", average internal = " << sum/((double) v.size()) << std::endl;
      t.next("stats");
    }
  }
  t.next("delete tree");
}
