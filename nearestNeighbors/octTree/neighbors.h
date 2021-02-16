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
// 0=root based, 2=bit based, >3=map based
// 1 is cursed
int queue_cutoff = 10; 


#include <algorithm>
#include <math.h> 
#include <queue>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "oct_tree.h"


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

  static bool compare(vtx_dist a, vtx_dist b){
    return a.second < b.second; 
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

     
    std::priority_queue<vtx_dist, std::vector<vtx_dist>, decltype(&(k_nearest_neighbors::compare))> nearest_nbh;      


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
      }
      max_distance = numeric_limits<double>::max();
      nearest_nbh.push(std::make_pair(p, 0));
      nearest_nbh.pop();
      std::cout << "made kNN" << std::endl; 
    }
    

    // if p is closer than neighbors[0] then swap it in
    void update_nearest(vtx *other) {  
      auto dist = (vertex->pt - other->pt).sqLength();
      if (dist < max_distance) { 
      	neighbors[0] = other;
      	distances[0] = dist;
        max_distance = dist; 
      	for (int i = 1;
      	     i < k && distances[i-1] < distances[i];
      	     i++) {
      	  swap(distances[i-1], distances[i]);
      	  swap(neighbors[i-1], neighbors[i]); 
        }
      }
    }

    //put into queue if vtx is closer than the furthest neighbor
    void update_nearest_queue(vtx* other){
      std::cout << "updating queue" << std::endl;
      // std::cout << (vertex->pt - other->pt).sqLength();
      auto dist = (vertex->pt - other->pt).sqLength();
      if (dist < max_distance){
        // std::cout << "popping" << std::endl; 
        nearest_nbh.pop();
        // std::cout << "popped, now pushing new elt" << std::endl; 
        nearest_nbh.push(std::make_pair(other, dist));
        // std::cout << "pushed" << std::endl; 
        max_distance = dist; 
      }
      std::cout << "successfully updated queue" << std::endl; 
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

    //this is not the smartest, but it is unlikely to be a bottleneck
    void merge_queue(kNN &L, kNN &R){
      while (nearest_nbh.size()>0) nearest_nbh.pop();
      int L_size = L.nearest_nbh.size();
      int R_size = R.nearest_nbh.size();
      int total = L_size + R_size;
      max_distance = max(L.max_distance, R.max_distance);   
      parlay::sequence<vtx_dist> queue_merge;
      queue_merge = parlay::sequence<vtx_dist>();
      for(int i=0; i<L_size; i++){
        queue_merge[i] = L.nearest_nbh.top();
        L.nearest_nbh.pop();
      }
      for(int i=0; i<total; i++){
        queue_merge[L_size+i] = R.nearest_nbh.top();
        R.nearest_nbh.pop();
      }
      auto compare = [&] (vtx_dist a, vtx_dist b){
            return a.second < b.second;
      };
      auto x = parlay::sort(queue_merge, compare);
      int i=0;
      vtx_dist prev;
      while (i<k){
        if(x[i] != prev){
          nearest_nbh.push(x[i]);
          prev = x[i];
          i += 1; 
        }
      }
    }
    
    // looks for nearest neighbors for pt in Tree node T
    void k_nearest_rec(node* T) {
      if (report_stats) internal_cnt++;
      if (within_epsilon_box(T, sqrt(max_distance))) { 
         // std::cout << "here4" << std::endl; 
	       if (T->is_leaf()) {
	         if (report_stats) leaf_cnt++;
	         auto &Vtx = T->Vertices();
	         for (int i = 0; i < T->size(); i++)
	           if (Vtx[i] != vertex){ 
                if (k < queue_cutoff){
                  update_nearest(Vtx[i]);
                } else{
                  // std::cout << "here5" << std::endl; 
                  update_nearest_queue(Vtx[i]);
                  // std::cout << "here6" << std::endl; 
                }
              } 
	} else if (T->size() > 10000 && algorithm_version != 0) { 
	  auto L = *this; // make copies of the distances
	  auto R = *this; // so safe to call in parallel
	  parlay::par_do([&] () {L.k_nearest_rec(T->Left());},
			 [&] () {R.k_nearest_rec(T->Right());});
    if (k < queue_cutoff){
	     merge(L,R); // merge the results
    } else{
      merge_queue(L, R);
    }
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


//   //instantiates the leaf_sequence as a member of k_nearest_neighbors
//   using node_index = std::pair<size_t, node*>;
//   parlay::sequence<node_index> leaf_sequence; 


//   void populate_leaf_sequence0(node* T, box b, double Delta, int min, int max){
//     if(T -> is_leaf()){
//       size_t lower_left = o_tree::interleave_bits(b.first, b.first, Delta); 
//       size_t upper_right = o_tree::interleave_bits(b.second, b.first, Delta);
//       size_t smallest;
//       if(lower_left < upper_right) smallest = lower_left;
//       else smallest = upper_right;
//       node_index current = std::make_pair(smallest, T);
//       leaf_sequence[min] = current;
//     } else{
//       node* L = T->Left();
//       node* R = T->Right();
//       int L_children = L->num_leaf_children;
//       int R_children = R->num_leaf_children;
//       parlay::par_do_if(max > 1000,
//         [&] () {populate_leaf_sequence0(L, b, Delta, min, min+L_children-1);},
//         [&] () {populate_leaf_sequence0(R, b, Delta, min+L_children, max);}
//       );
//     }
//   }

//   void populate_leaf_sequence(node* R, box b, double Delta){
//     int children = R->num_leaf_children;
//     leaf_sequence = parlay::sequence<node_index>(children);
//     populate_leaf_sequence0(R, b, Delta, 0, children-1);
//     auto less = [] (node_index a, node_index b){
//       return a.first < b.first;
//     };
//     parlay::sort(leaf_sequence);
//   }


// //given a point, finds a leaf by binary searching through the leaf_sequence
// node* find_leaf_using_seq0(point p, box b, double Delta){ 
//   //find the interleaved version of the point
//   size_t searchInt = o_tree::interleave_bits(p, b.first, Delta); 
//   auto less = [&] (node_index a){
//     return a.first < searchInt;
//   };
//   size_t closest_index = parlay::internal::binary_search(leaf_sequence, less); 
//   node* closest = leaf_sequence[closest_index].second;
//   return closest; 
// }

// node* find_leaf_using_seq(point p, box b, double Delta){ //pointer to the root so it can interleave the bits
//   node* T = find_leaf_using_seq0(p, b, Delta);
//   return T;
// }


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

}; //this ends the k_nearest_neighbors structure




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
    knn_tree T(v);
    t.next("build tree");

    if (report_stats) 
      std::cout << "depth = " << T.tree->depth() << std::endl;

    // *******************
    if (algorithm_version == 0) { // this is for starting from root 
      // this reorders the vertices for locality
      parlay::sequence<vtx*> vr = T.vertices();
      t.next("flatten tree");

      // find nearest k neighbors for each point
      parlay::parallel_for (0, 2, [&] (size_t i) {
	       T.k_nearest(vr[i], k);
         std::cout << "found neighbor for point " << i << std::endl; 
      }, 1);

      // for (size_t i=0; i<v.size(); i++){
      //   T.k_nearest(vr[i], k);
      // }

    // *******************
    // } else if (algorithm_version == 1) { // using find_leaf
    //   // this reorders the vertices for locality
    //   parlay::sequence<vtx*> vr = T.vertices();
    //   t.next("flatten tree");

    //   int dims = (v[0]->pt).dimension();  
    //   node* root = T.tree.get(); 
    //   size_t size = v.size();
    //   box_delta bd = T.get_box_delta(root, dims);
    //   //construct the sorted sequence which indexes leaves with their center points
    //   T.populate_leaf_sequence(root, bd.first, bd.second);
    //   t.next("populate leaf sequence");
    //   node* leaf = T.find_leaf_using_seq(vr[0]->pt, bd.first, bd.second);
    //   bool check = leaf->is_leaf(); 
    //   if(check == true) std::cout << "yes" << "\n";
    //   t.next("find leaf for last point");
      
    //   parlay::parallel_for(size-1, size, [&] (size_t i) {   
	   //     T.k_nearest_leaf(vr[i], T.find_leaf_using_seq(vr[i]->pt, bd.first, bd.second), k);
    //    }
    //    );

    // *******************      
    } else if (algorithm_version == 2) {
        parlay::sequence<vtx*> vr = T.vertices();
        t.next("flatten tree");

        int dims = (v[0]->pt).dimension();  
        node* root = T.tree.get(); 
        box_delta bd = T.get_box_delta(root, dims);
        size_t size = v.size();

        parlay::parallel_for(0, size, [&] (size_t i) {
          T.k_nearest_leaf(vr[i], T.find_leaf(vr[i]->pt, root, bd.first, bd.second), k);
        }
        );


    } else { //(algorithm_version == 3) this is for starting from leaf, finding leaf using map()
        auto f = [&] (vtx* p, node* n){ 
  	     return T.k_nearest_leaf(p, n, k); 
        };

        // find nearest k neighbors for each point
        T.tree -> map(f);
    }

    t.next("try all");
    if (report_stats) {
      auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
      size_t i = parlay::max_element(s) - s.begin();
      size_t sum = parlay::reduce(s);
      std::cout << "max internal = " << s[i] 
		<< ", average internal = " << sum/((double) v.size()) << std::endl;
      t.next("stats");
    }
  t.next("delete tree");


};
}

