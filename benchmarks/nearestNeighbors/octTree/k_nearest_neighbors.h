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

int queue_cutoff = 50;   



#include <algorithm>
#include <math.h> 
#include <queue>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "oct_tree.h"
#include "qknn.hpp"


// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vtx, int max_k>
struct k_nearest_neighbors {
  using vtx_dist = std::pair<vtx*, double>;
  using point = typename vtx::pointT;
  using fvect = typename point::vector;
  using o_tree = oct_tree<vtx>;
  using node = typename o_tree::node;
  using tree_ptr = typename o_tree::tree_ptr;
  using box = typename o_tree::box;
  using slice_t = typename o_tree::slice_t;

  tree_ptr tree;

  box tree_box; 

  void set_box(box b){
    tree_box = b; 
  }

    // generates the search structure
  k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    tree = o_tree::build(V); 
    node* root = tree.get();
    set_box(root->Box());
  }

  k_nearest_neighbors(parlay::sequence<vtx*> &V, box b) {
    //TODO add safety check
    box points_box = o_tree::get_box(V);
    int dims = V[0]->pt.dimension();
    bool ll_bad = false;
    bool ur_bad = false;
    for(int i=0; i<dims; i++){
      if(points_box.first[i] < b.first[i]){
        ll_bad = true;
      }
      if(points_box.second[i] > b.second[i]){
        ur_bad = true;
      }
    }
    if(ll_bad || ur_bad){
      std::cout << "ERROR: user-specified box does not contain dataset" << std::endl;
      abort();
    }
    set_box(b);
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
    double max_distance; // needed since we may need to update our biggest boi without a vector
    int k;
    int dimensions;
    size_t leaf_cnt;
    size_t internal_cnt;
    qknn<vtx> nearest_nbh;      


    kNN() {}


    // returns the ith smallest element (0 is smallest) up to k-1
    // no need to make a queue equivalent
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
      if (k < queue_cutoff){
        for (int i=0; i<k; i++) {
        	neighbors[i] = (vtx*) NULL; 
        	distances[i] = numeric_limits<double>::max();
        }
      } else{
        nearest_nbh = qknn<vtx>();
        nearest_nbh.set_size(k);
      }
      max_distance = numeric_limits<double>::max();
    }
    

    // if p is closer than neighbors[0] then swap it in
    void update_nearest(vtx *other) {  
      auto dist = (vertex->pt - other->pt).sqLength();
      if (dist < max_distance) { 
      	neighbors[0] = other;
      	distances[0] = dist;
      	for (int i = 1;
      	     i < k && distances[i-1] < distances[i];
      	     i++) {
      	  swap(distances[i-1], distances[i]);
      	  swap(neighbors[i-1], neighbors[i]); 
        }
        max_distance = distances[0];
      }
    }

    //put into queue if vtx is closer than the furthest neighbor
    void update_nearest_queue(vtx* other){
      auto dist = (vertex->pt - other->pt).sqLength();
      bool updated = nearest_nbh.update(other, dist);
      if (updated){
        max_distance = nearest_nbh.topdist();
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
      return (T->center() - vertex->pt).sqLength();
    }

   

    // sorted backwards
    void merge(kNN &L, kNN &R) {
      int i = k-1;
      int j = k-1;
      int r = k-1;
      while (r >= 0) {
	if (L.distances[i] < R.distances[j]) {
	  distances[r] = L.distances[i];
	  neighbors[r] = L.neighbors[i];
	  i--; 
	} else {
	  distances[r] = R.distances[j];
	  neighbors[r] = R.neighbors[j];
	  // same neighbor could appear in both lists
	  if (L.neighbors[i] == R.neighbors[j]) i--;
	  j--;
	}
	r--;
      }
    }
    
    // looks for nearest neighbors for pt in Tree node T
    void k_nearest_rec(node* T) {
      if (report_stats) internal_cnt++;
      if (within_epsilon_box(T, sqrt(max_distance))) { 
	       if (T->is_leaf()) {
	         if (report_stats) leaf_cnt++;
	         auto &Vtx = T->Vertices();
	         for (int i = 0; i < T->size(); i++)
	           if (Vtx[i] != vertex){ 
                if (k < queue_cutoff){
                  update_nearest(Vtx[i]);
                } else{
                  update_nearest_queue(Vtx[i]);
                }
              } 
	} else if (T->size() > 10000 && algorithm_version != 0 && k < queue_cutoff) { 
	  auto L = *this; // make copies of the distances
	  auto R = *this; // so safe to call in parallel
	  parlay::par_do([&] () {L.k_nearest_rec(T->Left());},
			 [&] () {R.k_nearest_rec(T->Right());});
	  merge(L,R); // merge the results
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
    if (current -> is_leaf()){
        if (report_stats) leaf_cnt++;
        auto &Vtx = T->Vertices();
        for (int i = 0; i < T->size(); i++)
          if (Vtx[i] != vertex){
            if (k < queue_cutoff){
              update_nearest(Vtx[i]);
            } else{
              update_nearest_queue(Vtx[i]);
            }
          }
      } 
    while((not within_epsilon_box(current, -sqrt(max_distance))) and (current -> Parent() != nullptr)){ 
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


  using box_delta = std::pair<box, double>;

  box_delta get_box_delta(int dims){
    box b = tree_box; 
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]);
    box_delta bd = make_pair(b, Delta);
    return bd;
  }


  // takes in an integer and a position in said integer and returns whether the bit at that position is 0 or 1
  int lookup_bit(size_t interleave_integer, int pos){ //pos must be less than key_bits, can I throw error if not?
    size_t val = ((size_t) 1) << (pos - 1);
    size_t mask = (pos == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << pos);
    if ((interleave_integer & mask) <= val){
      return 1;
    } else{
      return 0;
    }
  }

//This finds the leaf in the search structure that p is located in
node* find_leaf(point p, node* T, box b, double Delta){ //takes in a point since interleave_bits() takes in a point
  //first, we use code copied over from oct_tree to go from a point to an interleave integer
  node* current = T;
  size_t searchInt = o_tree::interleave_bits(p, b.first, Delta); //calling interleave_bits from oct_tree
  //then, we use this interleave integer to find the correct leaf
  while (not (current->is_leaf())){
    if(lookup_bit(searchInt, current -> bit) == 0){ 
      current = current->Right(); 
    } else{
      current = current->Left();
    }
  };
  return current;
}

//this instantiates the knn search structure and then calls the function k_nearest_fromLeaf
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

 
  parlay::sequence<vtx*> z_sort(parlay::sequence<vtx*> v, box b, double Delta){ 
    using indexed_point = typename o_tree::indexed_point; 
    size_t n = v.size();
    parlay::sequence<indexed_point> points;
    points = parlay::sequence<indexed_point>(n);
    parlay::parallel_for(0, n, [&] (size_t i){
      size_t p1 = o_tree::interleave_bits(v[i]->pt, b.first, Delta);
      indexed_point i1 = std::make_pair(p1, v[i]);
      points[i] = i1; 
    });
    auto less = [&] (indexed_point a, indexed_point b){
      return a.first < b.first;
    };
    auto x = parlay::sort(points, less);
    parlay::sequence<vtx*> v3; 
    v3 = parlay::sequence<vtx*>(n);
    parlay::parallel_for(0, n, [&] (size_t i){
      v3[i] = x[i].second; 
    });
    return v3; 
  }

    using indexed_point = typename o_tree::indexed_point; 

  indexed_point get_point(node* T){
    if (T->is_leaf()){
      return (T->indexed_pts)[0];
    } else{
      return get_point(T->Left());
    }
  }

  void batch_insert0(slice_t idpts, node* T){
    if (idpts.size()==0) return;
    T->set_flag(true);
    if(T->is_leaf()){
      if(T->size() + idpts.size() < o_tree::node_cutoff || T->bit == 0){
        T->batch_update(idpts);
      } else{    
        o_tree::batch_split(idpts, T);
      }
    } else{
        T->set_size(T->size()+idpts.size());
        int new_bit = T->bit; 
        size_t val = ((size_t) 1) << (new_bit - 1);
        size_t mask = (new_bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << new_bit);
        auto less = [&] (indexed_point x) {
          return (x.first & mask) < val;
        };
        int cut_point = parlay::internal::binary_search(idpts, less);
        int child_bit = (T->Right())->bit;
        if(child_bit == (new_bit-1)){
          parlay::par_do_if(idpts.size() > 100,
            [&] () {batch_insert0(idpts.cut(0, cut_point), T->Left());},
            [&] () {batch_insert0(idpts.cut(cut_point, idpts.size()), T->Right());}
          );
        } else{
          indexed_point sample = get_point(T);
          int sample_pos = lookup_bit(sample.first, new_bit-1);
          if(sample_pos == 0){ //the points already in the tree are on the left
            o_tree::create_new(T, idpts.cut(0, cut_point), new_bit, true);
            batch_insert0(idpts.cut(cut_point, idpts.size()), T->Left()); 
          } else{
            o_tree::create_new(T, idpts.cut(cut_point, idpts.size()), new_bit, false);
            batch_insert0(idpts.cut(0, cut_point), T->Right());
          } 
        }
    }
  }

  void batch_insert(parlay::sequence<vtx*> v, node* R, box b, double Delta){
    size_t vsize = v.size();
    //make sure all the points are within the bounding box of the initial tree
    box points_box = o_tree::get_box(v);
    int dims = v[0]->pt.dimension();
    bool ll_bad = false;
    bool ur_bad = false;
    for(int i=0; i<dims; i++){
      if(points_box.first[i] < b.first[i]){
        ll_bad = true;
      }
      if(points_box.second[i] > b.second[i]){
        ur_bad = true;
      }
    }
    if(ll_bad || ur_bad){
      std::cout << "ERROR: points not contained in bounding box of data structure" << std::endl;
      abort();
    }
    parlay::sequence<indexed_point> idpts; 
    idpts = parlay::sequence<indexed_point>(vsize);
    auto points = parlay::delayed_seq<indexed_point>(vsize, [&] (size_t i) -> indexed_point {
  return std::make_pair(o_tree::interleave_bits(v[i]->pt, b.first, Delta), v[i]);
      });
    auto less = [] (indexed_point a, indexed_point b){
      return a.first < b.first; 
    };
    auto x = parlay::sort(points, less);
    batch_insert0(parlay::make_slice(x), R);
    box root_box = o_tree::update_boxes(R);
  }

  void batch_delete0(slice_t idpts, node* R){
    if(idpts.size()==0) return;
    R->set_flag(true);
    if(R->is_leaf()){
      size_t n = idpts.size();
      if(n == R->size()){
        o_tree::prune(R);
      }else{
        R->batch_remove(idpts);
      }
    } else{
      size_t n = idpts.size();
      if(n == R->size()){
        o_tree::prune(R);
      } else{
        R->set_size(R->size()-n);
        int new_bit = R->bit; 
        size_t val = ((size_t) 1) << (new_bit - 1);
        size_t mask = (new_bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << new_bit);
        auto less = [&] (indexed_point x) {
          return (x.first & mask) < val;
        };
        int cut_point = parlay::internal::binary_search(idpts, less);
        parlay::par_do_if(n > 100,
          [&] () {batch_delete0(idpts.cut(0, cut_point), R->Left());},
          [&] () {batch_delete0(idpts.cut(cut_point, n), R->Right());}
        );
      }
    }
  }

  void batch_delete(parlay::sequence<vtx*> v, node* R, box b, double Delta){
    size_t vsize = v.size();
    parlay::sequence<indexed_point> idpts; 
    idpts = parlay::sequence<indexed_point>(vsize);
    auto points = parlay::delayed_seq<indexed_point>(vsize, [&] (size_t i) -> indexed_point {
      return std::make_pair(o_tree::interleave_bits(v[i]->pt, b.first, Delta), v[i]);
    });
    auto less = [] (indexed_point a, indexed_point b){
      return a.first < b.first; 
    };
    auto x = parlay::sort(points, less); 
    batch_delete0(parlay::make_slice(x), R);
    o_tree::compress(R);  
    box root_box = o_tree::update_boxes(R);
  }

}; //this ends the k_nearest_neighbors structure

