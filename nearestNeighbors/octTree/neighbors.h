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

bool report_stats = false;
int algorithm_version = 2;

#include <algorithm>
<<<<<<< HEAD
#include <math.h> // so we can have the square root
=======
#include <math.h> 
>>>>>>> Magdalen
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

  tree_ptr tree;

    // generates the search structure
  k_nearest_neighbors(parlay::sequence<vtx*> &V) {
<<<<<<< HEAD
    tree = o_tree::build(V); //double colon because o_tree is a static method, not specific to the tree
=======
    tree = o_tree::build(V); 
>>>>>>> Magdalen
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
<<<<<<< HEAD
	} else if (T->size() > 10000) { 
=======
	} else if (T->size() > 10000 && algorithm_version != 0) { 
>>>>>>> Magdalen
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


<<<<<<< HEAD
=======
  using box_delta = std::pair<box, double>;

  box_delta get_box_delta(node* T, int dims){
    box b = T -> Box(); 
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]);
    box_delta bd = make_pair(b, Delta);
    return bd;
  }


>>>>>>> Magdalen
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
<<<<<<< HEAD
node* find_leaf(point p, int dims, node* T){ //takes in a point since interleave_bits() takes in a point
  //first, we use code copied over from oct_tree to go from a point to an interleave integer
  using box = typename o_tree::box;
  node* current = T;
  box b = current -> Box(); 
  double Delta = 0;
  for (int i = 0; i < dims; i++) 
    Delta = std::max(Delta, b.second[i] - b.first[i]);
=======
node* find_leaf0(point p, node* T, box b, double Delta){ //takes in a point since interleave_bits() takes in a point
  //first, we use code copied over from oct_tree to go from a point to an interleave integer
  node* current = T;
>>>>>>> Magdalen
  size_t searchInt = o_tree::interleave_bits(p, b.first, Delta); //calling interleave_bits from oct_tree
  //then, we use this interleave integer to find the correct leaf
  while (not (current->is_leaf())){
    if(lookup_bit(searchInt, current -> bit) == 0){ 
      current = current->Right(); 
    } else{
      current = current->Left();
    }
<<<<<<< HEAD
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
}; 


//this instantiates the knn search structure and then calls the function k_nearest_fromLeaf
void k_nearest_leaf(vtx* p, node* T, int k) { 
  kNN nn(p, k); 
  nn.k_nearest_fromLeaf(T);
  if (report_stats) p->counter = nn.internal_cnt;
  for (int i=0; i < k; i++)
    p->ngh[i] = nn[i];
}

=======
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

node* find_leaf(point p, int dims, node* R){
  box_delta bd = get_box_delta(R, dims);
  node* T = find_leaf0(p, R, bd.first, bd.second);
  return T; 
}

  //instantiates the leaf_sequence as a member of k_nearest_neighbors
  using node_index = std::pair<size_t, node*>;
  parlay::sequence<node_index> leaf_sequence = parlay::sequence<node_index>(); 


  // sort the sequence based on the interlaved integers so we can binary search it later
  void sort_leaf_sequence(){
    auto less = [] (node_index a, node_index b) { 
      return a.first < b.first;
    };                                            
    parlay::sort(leaf_sequence, less); //TODO does using the center point always work?         
  }


  // for a node, add a pair of an interleaved center point and a node pointer to the leaf sequence
  void g_internal(node* T, box b, double Delta){
    size_t interleaved_center = o_tree::interleave_bits(T -> center(), b.first, Delta); 
    node_index current = std::make_pair(interleaved_center, T);
    leaf_sequence.push_back(current); //TODO check this is right
  }


  // populate the sequence with pairs of interleaved center points and node pointers
  void get_leaf_sequence0(node* R, box b, double Delta, size_t size){
    leaf_sequence.resize(size); //re-assign to just a sequence of the correct size
    // wrapper function so that map_node() can be called on g_internal
    auto g = [&] (node* T){ 
      return g_internal(T, b, Delta); 
    };
    R -> map_node(g);
  }

  void get_leaf_sequence(node* R, int dims, size_t size){
    box_delta bd = get_box_delta(R, dims);
    get_leaf_sequence0(R, bd.first, bd.second, size);
    sort_leaf_sequence();
  }

//given a point, finds a leaf by binary searching through the leaf_sequence
node* find_leaf_using_seq0(point p, box b, double Delta){ 
  //find the interleaved version of the point
  size_t searchInt = o_tree::interleave_bits(p, b.first, Delta); //TODO does center have only integer coords?
  auto less = [] (node_index a, size_t b){
    return a.first < b;
  };
  node_index closest = parlay::internal::binary_search(leaf_sequence, searchInt, less);
  return closest.second;
}

node* find_leaf_using_seq(point p, node* R, int dims){ //pointer to the root so it can interleave the bits
  box_delta bd = get_box_delta(R, dims); //TODO move this to ANN()
  node* T = find_leaf_using_seq0(p, bd.first, bd.second);
  return T;
}

void delete_leaf_seq(){ //TODO this is messed up
  leaf_sequence.clear();
}

//this instantiates the knn search structure and then calls the function k_nearest_fromLeaf
void k_nearest_leaf(vtx* p, node* T, int k) { 
  kNN nn(p, k); 
  nn.k_nearest_fromLeaf(T);
  if (report_stats) p->counter = nn.internal_cnt;
  for (int i=0; i < k; i++)
    p->ngh[i] = nn[i];
}

>>>>>>> Magdalen
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

<<<<<<< HEAD
    int type = 2;

    // *******************
    if (type == 0) { // this is for starting from root 
=======
    // *******************
    if (algorithm_version == 0) { // this is for starting from root 
>>>>>>> Magdalen
      // this reorders the vertices for locality
      parlay::sequence<vtx*> vr = T.vertices();
      t.next("flatten tree");

      // find nearest k neighbors for each point
      parlay::parallel_for (0, v.size(), [&] (size_t i) {
	  T.k_nearest(vr[i], k);}, 1);

    // *******************
<<<<<<< HEAD
    } else if (type == 1) { // using find_leaf
      // this reorders the vertices for locality
      parlay::sequence<vtx*> vr = T.vertices();
      t.next("flatten tree");
      
      int dims = (v[0]->pt).dimension();          
      parlay::parallel_for(0, v.size(), [&] (size_t i) {   
	  T.k_nearest_leaf(vr[i], T.find_leaf(vr[i]->pt, dims, T.tree.get()), k);});

    // *******************      
    } else { //(type == 2) this is for starting from leaf, finding leaf using map()
      auto f = [&] (vtx* p, node* n){ 
	return T.k_nearest_leaf(p, n, k); 
      };

      // find nearest k neighbors for each point
      T.tree -> map(f);
=======
    } else if (algorithm_version == 1) { // using find_leaf
      // this reorders the vertices for locality
      parlay::sequence<vtx*> vr = T.vertices();
      t.next("flatten tree");

      int dims = (v[0]->pt).dimension();   
      size_t size = v.size();
      //construct the sequence which indexes leaves with their center points, and sorts it
      T.get_leaf_sequence(T.tree.get(), dims, size);
      T.delete_leaf_seq();

      
             
      parlay::parallel_for(0, size, [&] (size_t i) {   
	  T.k_nearest_leaf(vr[i], T.find_leaf(vr[i]->pt, dims, T.tree.get()), k);});

    // *******************      
    } else { //(algorithm_version == 2) this is for starting from leaf, finding leaf using map()
        auto f = [&] (vtx* p, node* n){ 
  	     return T.k_nearest_leaf(p, n, k); 
        };

        // find nearest k neighbors for each point
        T.tree -> map(f);
>>>>>>> Magdalen
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

