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
#include "oct_tree.h"

// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vtx, int maxK>
struct k_nearest_neighbor {
  using point = typename vtx::pointT;
  using fvect = typename point::vector;

  using o_tree = oct_tree<vtx>;
  o_tree *tree;
  std::atomic<size_t> leaf_cnt = 0;
  std::atomic<size_t> internal_cnt = 0;
  
  // generates the search structure
  k_nearest_neighbor(parlay::sequence<vtx*> &V) {
    tree = o_tree::build(V);
  }

  // returns the vertices in the search structure, in an
  //  order that has spacial locality
  parlay::sequence<vtx*> vertices() {
    return tree->flatten();
  }

  struct kNN {
    vtx *ps;  // the element for which we are trying to find a NN
    vtx *pn[maxK];  // the current k nearest neighbors (nearest last)
    double rn[maxK]; // radius of current k nearest neighbors
    int k;
    int dim;
    size_t leaf_cnt;
    size_t internal_cnt;
    kNN() {}

    // returns the ith smallest element (0 is smallest) up to k-1
    vtx* operator[] (const int i) { return pn[k-i-1]; }

    kNN(vtx *p, int kk) {
      if (kk > maxK) {
	std::cout << "k too large in kNN" << std::endl;
	abort();}
      k = kk;
      ps = p;
      dim = p->pt.dimension();
      leaf_cnt = internal_cnt = 0;
      for (int i=0; i<k; i++) {
	pn[i] = (vtx*) NULL; 
	rn[i] = numeric_limits<double>::max();
      }
    }

    // if p is closer than pn then swap it in
    void update(vtx *p) { 
      point opt = (p->pt);
      fvect v = (ps->pt) - opt;
      double r = v.Length();
      if (r < rn[0]) {
	pn[0]=p; rn[0] = r;
	for (int i=1; i < k && rn[i-1]<rn[i]; i++) {
	  swap(rn[i-1],rn[i]); swap(pn[i-1],pn[i]); }
      }
    }

    bool within_epsilon_box(o_tree *T, double epsilon) {
      auto b = T->Box();
      bool result = true;
      for (int i = 0; i < dim; i++) {
	result = (result &&
		  (b.first[i] - epsilon < ps->pt[i]) &&
		  (b.second[i] + epsilon > ps->pt[i]));
      }
      return result;
    }

    double dist(o_tree *T) {
      return (T->center() - ps->pt).Length();
    }
    
    // looks for nearest neighbors for pt in Tree node T
    void nearestNgh(o_tree *T) {
      if (within_epsilon_box(T, rn[0])) {
	if (T->is_leaf()) {
	  leaf_cnt++;
	  for (int i = 0; i < T->size(); i++)
	  if (T->P[i] != ps) update(T->P[i]);
	} else if (dist(T->Left()) < dist(T->Right())) {
	  internal_cnt++;
	  nearestNgh(T->Left());
	  nearestNgh(T->Right());
	} else {
	  internal_cnt++;
	  nearestNgh(T->Right());
	  nearestNgh(T->Left());
	}
      }
    }
  };

  void k_nearest(vtx *p, int k) {
    kNN nn(p,k);
    nn.nearestNgh(tree);
    leaf_cnt += nn.leaf_cnt;
    internal_cnt += nn.internal_cnt;
    for (int i=0; i < k; i++)
      p->ngh[i] = nn[i];
  }

};

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN",false);
  using knn_tree = k_nearest_neighbor<vtx, maxK>;
  
  knn_tree T(v);
  t.next("build tree");
  size_t d = T.tree->depth();
  //std::cout << "depth = " << d << std::endl;

  // this reorders the vertices for locality
  parlay::sequence<vtx*> vr = T.vertices();
  t.next("flatten tree");
  
  // find nearest k neighbors for each point
  parlay::parallel_for (0, v.size(), [&] (size_t i) {
      T.k_nearest(vr[i], k);
    });
  //std::cout << "leaves = " << T.leaf_cnt.load()/((double) v.size()) 
  //          << ", internal = " << T.internal_cnt.load()/((double) v.size()) << std::endl;
  t.next("try all");
}
