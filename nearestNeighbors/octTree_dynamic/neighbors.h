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

bool report_stats = true;
int algorithm_version = 0;

#include <algorithm>
#include <math.h> 
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
  using box = typename o_tree::box;
  using slice_t = typename o_tree::slice_t;

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
      auto dist = (vertex->pt - other->pt).sqLength();
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
      if (within_epsilon_box(T, sqrt(distances[0]))) {
	if (T->is_leaf()) {
	  if (report_stats) leaf_cnt++;
	  auto &Vtx = T->Vertices();
	  for (int i = 0; i < T->size(); i++)
	    if (Vtx[i] != vertex) update_nearest(Vtx[i]);
	} else if (T->size() > 10000 && algorithm_version != 0) { 
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
          if (Vtx[i] != vertex) update_nearest(Vtx[i]);
      } 
    while((not within_epsilon_box(current, -sqrt(distances[0]))) and (current -> Parent() != nullptr)){ //check that current parent is not null
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

  box_delta get_box_delta(node* T, int dims){
    box b = T -> Box(); 
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
  //this is a test case, only works for a two-dimensional sample size
  // bool check = false; 
  // auto &Vtx = current -> Vertices();
  // for (int i = 0; i < current -> size(); i++)
  //   if ( ((Vtx[i] -> pt).x == p.x) and ((Vtx[i] -> pt).y == p.y)) { //this is a hack since it doesn't work properly for 3d
  //   //need a notion of equality here---check that components are equal?
  //     check = true;
  //   }
  // std::cout << "check " << check << "\n"; 
  return current;
}


  //instantiates the leaf_sequence as a member of k_nearest_neighbors
  using node_index = std::pair<size_t, node*>;
  parlay::sequence<node_index> leaf_sequence; 


  void populate_leaf_sequence0(node* T, box b, double Delta, int min, int max){
    if(T -> is_leaf()){
      size_t lower_left = o_tree::interleave_bits(b.first, b.first, Delta); 
      size_t upper_right = o_tree::interleave_bits(b.second, b.first, Delta);
      size_t smallest;
      if(lower_left < upper_right) smallest = lower_left;
      else smallest = upper_right;
      node_index current = std::make_pair(smallest, T);
      leaf_sequence[min] = current;
    } else{
      node* L = T->Left();
      node* R = T->Right();
      int L_children = L->num_leaf_children;
      int R_children = R->num_leaf_children;
      parlay::par_do_if(max > 1000,
        [&] () {populate_leaf_sequence0(L, b, Delta, min, min+L_children-1);},
        [&] () {populate_leaf_sequence0(R, b, Delta, min+L_children, max);}
      );
    }
  }

  void populate_leaf_sequence(node* R, box b, double Delta){
    int children = R->num_leaf_children;
    leaf_sequence = parlay::sequence<node_index>(children);
    populate_leaf_sequence0(R, b, Delta, 0, children-1);
    auto less = [] (node_index a, node_index b){
      return a.first < b.first;
    };
    parlay::sort(leaf_sequence);
  }


//given a point, finds a leaf by binary searching through the leaf_sequence
node* find_leaf_using_seq0(point p, box b, double Delta){ 
  //find the interleaved version of the point
  size_t searchInt = o_tree::interleave_bits(p, b.first, Delta); 
  auto less = [&] (node_index a){
    return a.first < searchInt;
  };
  size_t closest_index = parlay::internal::binary_search(leaf_sequence, less); 
  node* closest = leaf_sequence[closest_index].second;
  return closest; 
}

node* find_leaf_using_seq(point p, box b, double Delta){ //pointer to the root so it can interleave the bits
  node* T = find_leaf_using_seq0(p, b, Delta);
  return T;
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


  using indexed_point = typename o_tree::indexed_point; 

  indexed_point get_point(node* T){
    if (T->is_leaf()){
      return (T->indexed_pts)[0];
    } else{
      return get_point(T->Left());
    }
  }

  void batch_insert0(slice_t idpts, node* T){
    // std::cout << idpts.size() << std::endl; 
    if (idpts.size()==0) return;
    if(T-> is_leaf()){
      if(T->size() + idpts.size() < o_tree::node_cutoff || T->bit == 0){
        // std::cout << "updating node" << std::endl; 
        T->batch_update(idpts);
        // std::cout << "successfully updated" << std::endl; 
      } else{
        // std::cout << "splitting node" << std::endl;     
        o_tree::batch_split(idpts, T);
        // std::cout << "successfully split" << std::endl; 
      }
    } else{
        T->set_size(T->size()+idpts.size());
        // std::cout << "set new size" << std::endl; 
        int new_bit = T->bit; 
        size_t val = ((size_t) 1) << (new_bit - 1);
        size_t mask = (new_bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << new_bit);
        auto less = [&] (indexed_point x) {
          return (x.first & mask) < val;
        };
        int cut_point = parlay::internal::binary_search(idpts, less);
        int child_bit = (T->Right())->bit;
        // std::cout << "found cut point" << std::endl; 
        if(child_bit == (new_bit-1)){
          // std::cout << "recursing in normal case" << std::endl; 
          parlay::par_do_if(idpts.size() > 100,
            [&] () {batch_insert0(idpts.cut(0, cut_point), T->Left());},
            [&] () {batch_insert0(idpts.cut(cut_point, idpts.size()), T->Right());}
          );
        } else{
          indexed_point sample = get_point(T);
          int sample_pos = lookup_bit(sample.first, new_bit-1);
          // std::cout << "recursing in edge case" << std::endl; 
          if (sample_pos==0){
            o_tree::create_new(T, idpts.cut(0, cut_point), new_bit, true);
            batch_insert0(idpts.cut(cut_point, idpts.size()), T->Right());
          } else{
            o_tree::create_new(T, idpts.cut(cut_point, idpts.size()), new_bit, false);
            batch_insert0(idpts.cut(0, cut_point), T->Left());
          } 
        }
    }
  }

  void batch_insert(parlay::sequence<vtx*> v, node* R, box b, double Delta){
    size_t vsize = v.size();
    //make sure the points are all within the bounding box of the initial data structure
    size_t corner1 = o_tree::interleave_bits(b.first, b.first, Delta);
    size_t corner2 = o_tree::interleave_bits(b.second, b.first, Delta);
    parlay::parallel_for (0, vsize, [&] (size_t i) {
      size_t interleaved_point = o_tree::interleave_bits(v[i]->pt, b.first, Delta);
      bool oneway = (interleaved_point >= corner1 && interleaved_point <= corner2);
      bool otherway = (interleaved_point >= corner2 && interleaved_point <= corner1);
      if(not (oneway || otherway)){
        std::cout << "error: tried to insert point outside bounding box" << std::endl;
        abort(); 
      }
    }, 1);
    // std::cout << "checked bounding boxes" << std::endl; 
    parlay::sequence<indexed_point> idpts; 
    idpts = parlay::sequence<indexed_point>(vsize);
    // std::cout << "initialized sequences" << std::endl; 
    auto points = parlay::delayed_seq<indexed_point>(vsize, [&] (size_t i) -> indexed_point {
  return std::make_pair(o_tree::interleave_bits(v[i]->pt, b.first, Delta), v[i]);
      });
    // std::cout << "populated sequence" << std::endl; 
    auto less = [] (indexed_point a, indexed_point b){
      return a.first < b.first; 
    };
    auto x = parlay::sort(points, less);
    // std::cout << "entering recursive step" << std::endl; 
    batch_insert0(parlay::make_slice(x), R);
  }

}; //this ends the k_nearest_neighbors structure


template<class vtx>
void print_seq(parlay::sequence<vtx*> v){
  for(size_t i=0; i<v.size(); i++) std::cout << v[i] << std::endl; 
}

//the parameter p means that p points will be added dynamically
int p = 5000000;  

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int max_k, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN",report_stats);

  {
    using knn_tree = k_nearest_neighbors<vtx, max_k>;
    using node = typename knn_tree::node;
    using box = typename knn_tree::box;
    using box_delta = std::pair<box, double>;
    size_t n = v.size();
    size_t v_2 = n/2;
    size_t v_1 = n/2;
    // std::cout << v_2 << std::endl; 
    parlay::sequence<vtx*> v1;
    v1 = parlay::sequence<vtx*>(v_1);
    parlay::sequence<vtx*> v2;
    v2 = parlay::sequence<vtx*>(v_2);
    // std::cout << "initialized sequences" << std::endl; 
    parlay::parallel_for (0, n, [&] (size_t i) {
      if(i < n/2) v1[i] = v[i];
      else v2[i] = v[i];
    }, 1);
    t.next("made query sequences");
    knn_tree T(v1);
    t.next("build tree");
    int dims = v[0]->pt.dimension();
    node* root = T.tree.get();
    box_delta bd = T.get_box_delta(root, dims);
    // std::cout << "entering batch insertion" << std::endl; 
    
    // std::cout << "at batch insertion" << std::endl; 
      T.batch_insert(v2, root, bd.first, bd.second);
    t.next("do batch insertion");

    // if (report_stats) 
    //   std::cout << "depth = " << T.tree->depth() << std::endl;

  //   // *******************
  //   if (algorithm_version == 0) { // this is for starting from root 
  //     // this reorders the vertices for locality
  //     parlay::sequence<vtx*> vr = T.vertices();    
  //     t.next("flatten tree");
  //     // find nearest k neighbors for each point
  //     parlay::parallel_for (0, vr.size(), [&] (size_t i) {
	 //  T.k_nearest(vr[i], k);}, 1);

  //   // *******************
  //   } else if (algorithm_version == 2) {
  //       parlay::sequence<vtx*> vr = T.vertices();
  //       t.next("flatten tree");

  //       int dims = (v[0]->pt).dimension();  
  //       node* root = T.tree.get(); 
  //       box_delta bd = T.get_box_delta(root, dims);
  //       size_t size = v.size();

  //       parlay::parallel_for(0, size, [&] (size_t i) {
  //         T.k_nearest_leaf(vr[i], T.find_leaf(vr[i]->pt, root, bd.first, bd.second), k);
  //       }
  //       );


  //   }

  //   // *******************      
  //    else { //(algorithm_version == 3) this is for starting from leaf, finding leaf using map()
  //       auto f = [&] (vtx* p, node* n){ 
  // 	     return T.k_nearest_leaf(p, n, k); 
  //       };

  //       // find nearest k neighbors for each point
  //       T.tree -> map(f);
  //   }

  //   t.next("try all");
  //   if (report_stats) {
  //     auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
  //     size_t i = parlay::max_element(s) - s.begin();
  //     size_t sum = parlay::reduce(s);
  //     std::cout << "max internal = " << s[i] 
		// << ", average internal = " << sum/((double) v.size()) << std::endl;
  //     t.next("stats");
  //   }
  t.next("delete tree");


};
}

