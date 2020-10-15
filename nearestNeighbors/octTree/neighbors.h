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

#define report_stats true
#include <iostream>
#include <algorithm>
#include <math.h> // so we can have the square root
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "oct_tree.h"


// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vtx, int max_k>
struct k_nearest_neighbors {
  using point = typename vtx::pointT;
  using fvect = typename point::vector;
  using o_tree = oct_tree<vtx>;
  using node = typename o_tree::node;
  using tree_ptr = typename o_tree::tree_ptr;

  tree_ptr tree;

  // generates the search structure
  k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    tree = o_tree::build(V); //double colon because o_tree is a static method, not specific to the tree
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




    
    // looks for nearest neighbors for pt in Tree node T
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

  void k_nearest_fromLeaf(node* T) {
    node* current = T; //this will be the node that node*T points to
    if (current -> is_leaf()) { //move this part outside the while loop
        if (report_stats) leaf_cnt++;
        auto &Vtx = T->Vertices();
        for (int i = 0; i < T->size(); i++)
          if (Vtx[i] != vertex) update_nearest(Vtx[i]);
      } 
    while(not within_epsilon_box(current, -sqrt(distances[0]))){
      node* parent = (current -> Parent());
      if (current == parent -> Right()){
        k_nearest_rec(parent -> Left());
      } else{
        k_nearest_rec(parent -> Right());
      }
      current = parent;
      
    }

  }

  }; // this ends the knn structure

  static size_t interleave_bits(point p, point min_point, double delta) {
    int dim = p.dimension();
    int bits = 64/dim; //because key_bits == 64
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

  int lookup_bit(size_t interleave_integer, int pos){ //pos must be less than key_bits, can I throw error if not?
    size_t val = ((size_t) 1) << (pos - 1);
    size_t mask = (pos == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << pos);
    if ((interleave_integer & mask) <= val){
      return 1;
    } else{
      return 0;
    };
  }

//This finds the leaf in the search structure that p is located in

node* find_leaf(point p, parlay::sequence<vtx*> V, node* T){ //takes in a point since interleave_bits() takes in a point
  //first, we use code copied over from oct_tree to go from a point to an interleave integer
  using box = typename o_tree::box;
  node* current = T;
  int dims = (V[0]->pt).dimension();
  box b = current -> Box(); //need to get the type for box
  double Delta = 0;
  for (int i = 0; i < dims; i++) 
    Delta = std::max(Delta, b.second[i] - b.first[i]);
  size_t searchInt = interleave_bits(p, b.first, Delta); //calling interleave_bits from oct_tree
  //then, we use this interleave integer to find the correct leaf
  while (!(current->is_leaf())){
    if(lookup_bit(searchInt, current -> bit) == 0){ 
      current = current->Left(); //
    } else{
      current = current->Right();
    }
  };
  return current;
}; 

//Box() is a node attribute
//interleave_bits and lookup_bit are oct_tree functions
//bit is a node attribute

//this instantiates the knn search structure and then calls the function k_nearest_fromLeaf
//it should be handed to map()

//since the kNN is created for each point, creating the kNN must be handed to map

void k_nearest_leaf(vtx* p, node* T, int k) { 
  kNN nn(p, k); 
  nn.k_nearest_fromLeaf(T);
  if (report_stats) p->counter = nn.internal_cnt;
  for (int i=0; i < k; i++)
    p->ngh[i] = nn[i];
}

  void k_nearest(vtx* p, int k) {
    kNN nn(p,k);
    nn.k_nearest_rec(tree.get()); //this is passing in a pointer to the o_tree root
    if (report_stats) p->counter = nn.internal_cnt;
    for (int i=0; i < k; i++)
      p->ngh[i] = nn[i];
  }

}; //this ends the k_nearest_neighbors structure




// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int max_k, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN",report_stats);

  {
    using knn_tree = k_nearest_neighbors<vtx, max_k>;
    using node = typename knn_tree::node;
    knn_tree T(v);
    t.next("build tree");

    if (report_stats) 
      std::cout << "depth = " << T.tree->depth() << std::endl;

    // // this reorders the vertices for locality
    // parlay::sequence<vtx*> vr = T.vertices();
    // t.next("flatten tree");

    //get a node pointer to the root instead of a tree_ptr to match specs of find_leaf()
    node* rightChild = T.tree -> Right();
    node* rootptr2 = rightChild -> Parent();

    parlay::parallel_for(0, v.size(), [&] (size_t i) {
      T.k_nearest_leaf(v[i], T.find_leaf(v[i]->pt, v, rootptr2), k);});

    // auto f = [&] (vtx* p, node* n){ 
    //   return T.k_nearest_leaf(p, n, k); //fix this to have access to the function k_nearest_leaf
    // };

  
    // // find nearest k neighbors for each point
    // T.tree -> map(f); 

    // // find nearest k neighbors for each point
    // parlay::parallel_for (0, v.size(), [&] (size_t i) {
    //        T.k_nearest(vr[i], k);}, 1);

  //   t.next("try all");
  //   if (report_stats) {
  //     auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
  //     size_t i = parlay::max_element(s) - s.begin();
  //     size_t sum = parlay::reduce(s);
  //     std::cout << "max internal = " << s[i] 
		// << ", average internal = " << sum/((double) v.size()) << std::endl;
  //     t.next("stats");
  //   }
  // t.next("delete tree");
};
}

