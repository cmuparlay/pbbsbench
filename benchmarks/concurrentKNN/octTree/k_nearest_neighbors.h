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
#include "flock/flock.h"
#include "verlib/verlib.h"
#include <assert.h>
#include <sstream>



// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vtx, int max_k>
struct k_nearest_neighbors {
  using vtx_dist = std::pair<vtx*, double>;
  using point = typename vtx::pointT;
  using fvect = typename point::vector;
  using o_tree = oct_tree<vtx>;
  using node = typename o_tree::node;
  // using tree_ptr = typename o_tree::tree_ptr;
  using box = typename o_tree::box;
  using slice_t = typename o_tree::slice_t;

  verlib::versioned_ptr<node> tree;
  flck::lock root_lock;

  box tree_box; 

  bool box_eq(box b, box c, int d){
    bool first = true;
    bool second = true;
    for(int i=0; i<d; i++){
      first = first && (b.first[i] == c.first[i]);
      second = second && (b.second[i] == c.second[i]);
    }
    return (first && second);
  }

	static double dist_sq_to_box(box b, point p){
    if(within_box(b, p)) return 0;
		double total=0;
    int d = p.dimension();
    for(int i=0; i<d; i++){
      auto d1 = b.first[i]-p[i];
      auto d2 = p[i]-b.second[i];
      auto dist = std::max(std::max(d1, 0.0), std::max(d2, 0.0));
      total += dist*dist;
    }
    return total;
	}

  void are_equal(node* T, int d){
    node* V = tree.load();
    return are_equal_rec(V, T, d);
  }

  void are_equal_rec(node* V, node* T, int d){
    if(T->bit != V->bit){
      std::cout << "UNEQUAL: bit" << std::endl;
      std::cout << "Inserted tree has bit " << V->bit << " while regular tree has bit " << T->bit << std::endl;
    }
    if(!box_eq(T->Box(), V->Box(), d)){
      std::cout << "UNEQUAL: box" << std::endl;
    }
    if(!(T->is_leaf()) && !(V->is_leaf())){
      are_equal_rec(T->Left(), V->Left(), d);
      are_equal_rec(T->Right(), V->Right(), d);
      return;
    }
    else if(T->is_leaf() && V->is_leaf()){
      if(T->size() != V->size()){
        std::cout << "UNEQUAL: leaf size" << std::endl;
        std::cout << "Inserted tree has leaf size " << V->size() << " while regular tree has leaf size " << T->size() << std::endl;
        std::cout << "Leaves have bit " << V->bit << std::endl;
      } //not a true eq check
      return;
    } else{
      std::cout << "UNEQUAL: internal node vs leaf node" << std::endl;
      abort();
    }
  }

  void set_box(box &b){tree_box = b;}

  void set_root(node* root){tree.store(root);}

  void set_sizes(){
    node* T = tree.load();
    set_sizes_rec(T);
  }

  size_t set_sizes_rec(node* T){
    if(T->is_leaf()) return T->size();
    else{
      size_t s = set_sizes_rec(T->Left()) + set_sizes_rec(T->Right());
      T->set_size(s);
      return s;
    }
  }

    // generates the search structure
  k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    tree = o_tree::build(V); 
    set_box(tree.load()->Box());
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
    tree = o_tree::build(V, b); 
  }

  ~k_nearest_neighbors() {
    node::delete_rec(tree.load());
  }


  // returns the vertices in the search structure, in an
  //  order that has spacial locality
  parlay::sequence<vtx*> vertices() {
    return tree.load()->flatten();
  }

  template<typename A, typename B, typename C>
  struct tuple {
    bool a; bool b; bool c;
    tuple(bool a, bool b, bool c) : a(a), b(b), c(c) {}
    tuple(size_t val) : a(val&(1<<2)), b(val&(1<<1)), c(val&(1)) {}
    // operator size_t&() { return 0+a; }
    operator size_t() const { return (((size_t)a)<<2) + (((size_t)b)<<1) + (size_t)c; }
  };

  tuple<bool,bool,bool> make_tuple(bool a, bool b, bool c) { return tuple<bool,bool,bool>(a,b,c); }

  template<typename A, typename B>
  struct pair {
    bool first; bool second;
    pair(bool a, bool b) : first(a), second(b) {}
    pair(size_t val) : first(val&(1<<1)), second(val&(1)) {}
    // operator size_t&() { return 0+first; }
    operator size_t() const { return (((size_t)first)<<1) + (size_t)second; }
  };

  pair<bool,bool> make_pair(bool a, bool b) { return pair<bool,bool>(a,b); }

  struct kNN {
    vtx *vertex;  // the vertex for which we are trying to find a NN
    vtx *neighbors[max_k];  // the current k nearest neighbors (nearest last)
    double distances[max_k]; // distance to current k nearest neighbors
    double max_distance; 
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
        abort();
      }
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

    bool is_deep(node* T){
      return (T->Left()->depth() >= 20) && (T->Right()->depth() >= 20);
    }
    
    // looks for nearest neighbors for pt in Tree node T
    void k_nearest_rec(node* T) {
      box b;
      b = T->Box();
      if (within_epsilon_box(T, sqrt(max_distance))) { 
        if (report_stats) internal_cnt++;
        if (T->is_leaf()) {
          if (report_stats) leaf_cnt++;
          auto &Vtx = T->Indexed_Pts();
          for (int i = 0; i < T->size(); i++){
            if (Vtx[i].second != vertex) { 
              if (k < queue_cutoff){
                update_nearest(Vtx[i].second);
              } else{
                update_nearest_queue(Vtx[i].second);
            }
          } 
          }
        } 
        else if (false && T->size() >= 10000 && max_distance != numeric_limits<double>::max()) { // skipping this case because par_do does not interact well with snapshots
          auto L = *this; // make copies of the distances
          auto R = *this; // so safe to call in parallel
          L.k_nearest_rec(T->Left());
          R.k_nearest_rec(T->Right());
	  //parlay::par_do([&] () {L.k_nearest_rec(T->Left());},
          //  [&] () {R.k_nearest_rec(T->Right());});
          merge(L,R); // merge the results
        }
        else if (dist_sq_to_box(T->Left()->Box(), vertex->pt) < dist_sq_to_box(T->Right()->Box(), vertex->pt)) {
          k_nearest_rec(T->Left());
          k_nearest_rec(T->Right());
        } else {
          k_nearest_rec(T->Right());
          k_nearest_rec(T->Left());
        }
      }
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

}; // this ends the knn structure

struct RNG {
    vtx *vertex;  // the vertex for which we are trying to find a NN
    parlay::sequence<int> range_candidates;
    double radius; 
    double r_sq;
    int dimensions;
    size_t leaf_cnt;
    size_t internal_cnt;


    RNG() {}


    RNG(vtx *p, double rad) {
      vertex = p;
      radius = rad;
      r_sq = radius*radius;
      dimensions = p->pt.dimension();
      leaf_cnt = internal_cnt = 0;
    }

    parlay::sequence<int> return_answer(){
      return range_candidates;
    }

    bool within_radius(node* T) {
      auto box = T->Box();
      bool result = true;
      for (int i = 0; i < dimensions; i++) {
        result = (result &&
              (box.first[i] - radius < vertex->pt[i]) &&
              (box.second[i] + radius > vertex->pt[i]));
      }
      return result;
    }

    bool point_within_radius(vtx* other){
      return ((vertex->pt-other->pt).sqLength() <= r_sq);
    }
    
    // looks for nearest neighbors for pt in Tree node T
    void range_search_rec(node* T) {
      box b;
      b = T->Box();
      if (within_radius(T)) { 
        if (report_stats) internal_cnt++;
        if (T->is_leaf()) {
          if (report_stats) leaf_cnt++;
          auto &Vtx = T->Indexed_Pts();
          for (int i = 0; i < T->size(); i++){
            if (Vtx[i].second != vertex) { 
              if (point_within_radius(Vtx[i].second)) {range_candidates.push_back(Vtx[i].second->identifier);}
            } 
          }
        } else {
          range_search_rec(T->Left());
          range_search_rec(T->Right());
        }
      }
    }

}; // end RNG


  using box_delta = std::pair<box, double>;

  box_delta get_box_delta(int dims){
    box b = tree_box; 
    double Delta = 0;
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]);
    box_delta bd = std::make_pair(b, Delta);
    return bd;
  }


  // takes in an integer and a position in said integer and returns whether the bit at that position is 0 or 1
  static int lookup_bit(size_t interleave_integer, int pos){ //pos must be less than key_bits, can I throw error if not?
    size_t val = ((size_t) 1) << (pos - 1);
    size_t mask = (pos == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << pos);
    if ((interleave_integer & mask) < val){
      return 0;
    } else{
      return 1;
    }
  }

//This finds the leaf in the search structure that p is located in
node* find_leaf(node* T){ //takes in a point since interleave_bits() takes in a point
  //first, we use code copied over from oct_tree to go from a point to an interleave integer
  node* current = T;
  //then, we use this interleave integer to find the correct leaf
  while (!(current->is_leaf())){
      current = current->Left(); 
  }
  std::cout << current->Indexed_Pts().size() << std::endl;
  std::cout << current->size() << std::endl;
  auto &Vtx = current->Indexed_Pts();
  for (int i = 0; i < current->size(); i++){
    std::cout << Vtx[i].second->identifier << std::endl;

  }
  return current;
}

  int k_nearest(vtx* p, int k) {
    int vv;
    verlib::with_snapshot([&] {
      kNN nn(p,k);
      nn.k_nearest_rec(tree.load()); 
      vv = nn.internal_cnt;
    });
    return vv;
  }

  std::tuple<int, int, parlay::sequence<int>> range_search(vtx* p, double rad){
      RNG rn(p,rad);
      int ic;
      int lc;
      verlib::with_snapshot([&]{
        rn.range_search_rec(tree.load());
        lc = rn.internal_cnt;
        ic = rn.leaf_cnt;
      });
      return {ic, lc, std::move(rn.return_answer())};
  }
  

 
  static parlay::sequence<vtx*> z_sort(parlay::sequence<vtx*> v, box b){ 
    using indexed_point = typename o_tree::indexed_point; 
    double Delta = 0;
    int dims = v[0]->pt.dimension();
    for (int i = 0; i < dims; i++) 
      Delta = std::max(Delta, b.second[i] - b.first[i]);
    size_t n = v.size();
    parlay::sequence<indexed_point> points;
    points = parlay::sequence<indexed_point>(n);
    parlay::parallel_for(0, n, [&] (size_t i){
      size_t p1 = o_tree::interleave_bits(v[i]->pt, b.first, Delta);
      points[i] = std::make_pair(p1, v[i]);
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

  static bool within_box(node* T, vtx* vertex) {
    return within_box(T->Box(), vertex->pt);
  }

  static bool within_box(box b, point p){
    int dimensions = p.dimension();
    bool result = true;
    for (int i = 0; i < dimensions; i++) {
      result = (result &&
        (b.first[i] <= p[i]) &&
        (b.second[i] >= p[i]));
    }
    return result;
  }

  //strips a point of its child pointers before deleting
  //to ensure child pointers are not deleted
  void delete_single(node* T){
    if(T != nullptr) {
      T->removed.store(true);
      node::retire_node(T);
    }
  }

  void delete_single_no_retire(node* T){
    if(T != nullptr) {
      T->removed.store(true);
      // node::retire_node(T);
    }
  }

  //takes in a single point and a pointer to the root of the tree
  //as well as a bounding box and its largest side length
  //assumes integer coordinates
  bool insert_point(vtx* p){
    return verlib::with_epoch([&] {
      int dims = p->pt.dimension();
      auto [b, Delta] = get_box_delta(dims);
      size_t interleave_int = o_tree::interleave_bits(p->pt, b.first, Delta);
      indexed_point q = std::make_pair(interleave_int, p);
      while(true) {
        node* R = tree.load();
        // set_root(insert_point_path_copy_helper(q, R, o_tree::key_bits));
        // return true;
#ifdef PathCopy
        pair<bool, bool> result = insert_point_path_copy(q, nullptr, R, o_tree::key_bits);
#else
        pair<bool, bool> result = insert_point0(q, nullptr, R, o_tree::key_bits);
#endif
        // assert(result.second);
        // assert(result.first);
        if(result.second) return result.first;
        // std::cout << "insert attempt failed" << std::endl; 
      }
      // while(!insert_point0(q, nullptr, R, o_tree::key_bits, false)) {} // repeat until success
    });
  }

  struct link {
    node* n;
    bool go_left;
  };

  pair<bool, bool> insert_point_path_copy(indexed_point q, node* parent, node* T, int bit) {
    // if(true) {
    if(T->is_leaf() || !within_box(T, q.second)) {
      // enter locking mode
      if(parent == nullptr) { // T is root
        auto result = root_lock.try_lock_result([=] { 
          if(tree.load() != T) return make_pair(false, false);
          else {
            std::vector<link> path; path.reserve(4);
            path.push_back((link){nullptr, false});
            bool success = insert_point_path_copy_helper(q, T, bit, path);
            return make_pair(success, true);
          }
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); 
      } else {
        auto result = parent->lck.try_lock_result([=] {
          if(parent->removed.load() || !(parent->Left() == T || parent->Right() == T)) { // if P is removed or P's child isn't T
            // std::cout << "failed valdiation" << std::endl; 
            return make_pair(false, false);
          } else {
            std::vector<link> path; path.reserve(4);
            bool left_child = (parent->Left() == T);
            path.push_back((link){parent, left_child});
            bool success = insert_point_path_copy_helper(q, T, bit, path);
            return make_pair(success, true);
          } 
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); 
      }
    } else {
      bool go_left = (lookup_bit(q.first, T->bit) == 0);
      return insert_point_path_copy(q, T, 
              (go_left ? T->Left() : T->Right()), T->bit-1);
    }   
  } 

  void install_new_path(std::vector<link> &old_path, node* new_n) {
    for(int i = old_path.size()-1; i>= 1; i--) {
      node* T = old_path[i].n;
      box b = T->Box();
      
      if(old_path[i].go_left) {
        box bigger = box(T->Right()->Box().first.minCoords(new_n->Box().first), 
                       T->Right()->Box().second.maxCoords(new_n->Box().second));
        new_n = node::new_node(new_n, T->Right(), T->bit, bigger);
      }
      else {
        box bigger = box(T->Left()->Box().first.minCoords(new_n->Box().first), 
                       T->Left()->Box().second.maxCoords(new_n->Box().second));
        new_n = node::new_node(T->Left(), new_n, T->bit, bigger);
      }
    }
    if(old_path[0].n == nullptr)
      set_root(new_n);
    else
      old_path[0].n->set_child(new_n, old_path[0].go_left);
    for(int i = old_path.size()-1; i>= 1; i--)
      delete_single(old_path[i].n);
  }

  void install_new_path_insert(std::vector<link> &old_path, node* new_n) {
    for(int i = old_path.size()-1; i>= 1; i--) {
      node* T = old_path[i].n;
      box b = T->Box();
      
      if(old_path[i].go_left) {
        box bigger = box(T->Right()->Box().first.minCoords(new_n->Box().first), 
                       T->Right()->Box().second.maxCoords(new_n->Box().second));
        new_n = node::new_node(new_n, T->Right(), T->bit, bigger);
      }
      else {
        box bigger = box(T->Left()->Box().first.minCoords(new_n->Box().first), 
                       T->Left()->Box().second.maxCoords(new_n->Box().second));
        new_n = node::new_node(T->Left(), new_n, T->bit, bigger);
      }
    }
    if(old_path[0].n == nullptr)
      set_root(new_n);
    else
      old_path[0].n->set_child(new_n, old_path[0].go_left);
    for(int i = old_path.size()-1; i>= 1; i--)
      delete_single(old_path[i].n);
  }

  bool insert_point_path_copy_helper(indexed_point q, node* T, int bit, std::vector<link> &path) {
    if(T->is_leaf()) {
      for(indexed_point p : T->Indexed_Pts()) {
        if(points_equal(p.second->pt, q.second->pt)) return false;
      }
      std::vector<indexed_point> points;
      for(auto i : T->Indexed_Pts()) points.push_back(i);
      points.push_back(q);
      // auto points = parlay::tabulate(T->size()+1, [&] (long i) {
      //   return (i == T->size()) ? q : T->Indexed_Pts()[i];
      // }, 1000);
      //two cases: either (1) the leaf size is below the cutoff, in which case
      //we create a new leaf with the point q added, or (2) the leaf needs to be
      //split into an internal node and two leaf children
      if(T->size() + 1 < o_tree::node_cutoff || T->bit == 0){
        //case 1
        box b = T->Box();
        box bigger = box(b.first.minCoords(q.second->pt), b.second.maxCoords(q.second->pt));
        node* new_l = node::new_leaf(std::move(points), T->bit, bigger);
        install_new_path_insert(path, new_l);
        delete_single(T);
        return true;
      } else {
        //sort points in leaf by interleave order
        int cut_point = 0;
        auto less_sort = [&] (indexed_point a, indexed_point b){
          return a.first < b.first;
        };
        std::sort(points.begin(), points.end(), less_sort);
        //search for cut point
        int new_bit = T->bit;
        while((cut_point == 0 | cut_point == points.size()) && new_bit != 0){
          size_t val = ((size_t) 1) << (new_bit - 1);
          size_t mask = (new_bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << new_bit);
          auto less = [&] (indexed_point x) {
            return (x.first & mask) < val;
          };
          cut_point = parlay::internal::binary_search(points, less);
          new_bit--;
        }
        std::vector<indexed_point> left_s;
        std::vector<indexed_point> right_s;
        for(size_t i=0; i<points.size(); i++){
          if(i<cut_point) left_s.push_back(points[i]);
          else right_s.push_back(points[i]);
        }
        node* L = node::new_leaf(std::move(left_s), new_bit);
        node* R = node::new_leaf(std::move(right_s), new_bit);
        node* new_l = node::new_node(L, R, new_bit+1);
        install_new_path_insert(path, new_l);
        delete_single(T);
        return true;
      }
    } else { // expand bounding box
      return T->lck.with_lock([=] {
        auto path_mutable = path;
        bool go_left = (lookup_bit(q.first, T->bit) == 0);
        //two cases: (1) need to create new leaf and internal node,
        //or (2) bit unchanged and box changed

        //if next bit of integer isn't same as node,
        //iterate until either make a leaf or bits match
        node* Q = T;
        while(!(Q->is_leaf())){ Q=Q->Right(); }
        indexed_point s = Q->Indexed_Pts()[0];
        int cur_bit = bit;
        // std::cout << cur_bit << " " << T->bit << endl;
        assert(cur_bit >= T->bit);
        // if(parent != nullptr) assert(cur_bit < G->bit);
        while(cur_bit > T->bit) {
          if(lookup_bit(q.first, cur_bit) != lookup_bit(s.first, cur_bit)){
            // std::cout << "case 1" << std::endl;
            //we know we are in case 1
            //form leaf
            std::vector<indexed_point> points = {q};
            node* R = node::new_leaf(std::move(points), cur_bit-1);
            //new parent node should replace T as G's child
            node* P;
            if(lookup_bit(q.first, cur_bit) == 0) P = node::new_node(R, T, cur_bit);
            else P = node::new_node(T, R, cur_bit);
            install_new_path_insert(path_mutable, P);
            return true;
          } else cur_bit--;
        }
        // std::cout << "case 2" << std::endl;
        //case 2
        //calculate new box around T, and create new internal node
        //with corrected box
        node* child = (go_left ? T->Left() : T->Right());
        link lnk = (link){T, go_left};
        path_mutable.push_back(lnk);
        return insert_point_path_copy_helper(q, child, T->bit-1, path_mutable);
      });
    }
  }

  // return value: <q_added, locks successful>
  pair<bool, bool> insert_point0(indexed_point q, node* parent, node* T, int bit, bool parent_locked=false){
    //lock T
    if(!parent_locked && (T->is_leaf() || !within_box(T, q.second))) {
      // enter locking mode
      if(parent == nullptr) { // T is root
        auto result = root_lock.try_lock_result([=] { 
          if(tree.load() != T) return make_pair(false, false);
          else return insert_point0(q, parent, T, bit, true); 
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); 
      } else {
        auto result = parent->lck.try_lock_result([=] {
          if(parent->removed.load() || !(parent->Left() == T || parent->Right() == T)) { // if P is removed or P's child isn't T
            // std::cout << "failed valdiation" << std::endl; 
            return make_pair(false, false);
          } else return insert_point0(q, parent, T, bit, true);
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); 
      }
    }

    if(T->is_leaf()) {
      return make_pair(insert_into_leaf(q, parent, T), true);
    } else {
      if(!within_box(T, q.second)) {
        assert(parent_locked);
        return T->lck.with_lock([=] {
          //two cases: (1) need to create new leaf and internal node,
          //or (2) bit unchanged and box changed

          //if next bit of integer isn't same as node,
          //iterate until either make a leaf or bits match
          node* Q = T;
          while(!(Q->is_leaf())){Q=Q->Right();}
          indexed_point s = Q->Indexed_Pts()[0];
          int cur_bit = bit;
          assert(cur_bit >= T->bit);
          // if(parent != nullptr) assert(cur_bit < G->bit);
          while(cur_bit > T->bit) {
            if(lookup_bit(q.first, cur_bit) != lookup_bit(s.first, cur_bit)){
              // std::cout << "case 1" << std::endl;
              //we know we are in case 1
              //form leaf
              node* G = parent;
              std::vector<indexed_point> points = {q};
              node* R = node::new_leaf(points, cur_bit-1);
              //new parent node should replace T as G's child
              node* P;
              if(lookup_bit(q.first, cur_bit) == 0) P = node::new_node(R, T, cur_bit);
              else P = node::new_node(T, R, cur_bit);
              if(G != nullptr){
                bool left = false;
                if(T == G->Left()) left = true;
                G->set_child(P, left);
              }
              if(T == tree.load()) {
                // std::cout << "root changed (internal: added new level)" << std::endl; 
                set_root(P);
              }
              return make_pair(true, true);
            } else cur_bit--;
          }
          // std::cout << "case 2" << std::endl;
          //case 2
          //calculate new box around T, and create new internal node
          //with corrected box
          node* P = parent;
          node* L = T->Left();
          node* R = T->Right();
          box b = T->Box();
          box bigger = box(b.first.minCoords(q.second->pt), b.second.maxCoords(q.second->pt));
          node* N = node::new_node(L, R, T->bit, bigger);
          return N->lck.with_lock([=] {
            assert(tree.load() == T || P != nullptr);
            if(P != nullptr){
              if(T==P->Left()) P->set_child(N, true);
              else P->set_child(N, false);
#ifdef HandOverHand
              P->lck.unlock();
#endif
            } 
            if(T == tree.load()) {
              // std::cout << "root changed (internal: updated bounding box)" << std::endl;
              assert(parent == nullptr); 
              set_root(N);
#ifdef HandOverHand
              root_lock.unlock();
#endif
            }
            delete_single(T);
#ifdef HandOverHand
            T->lck.unlock();
#endif
            return insert_internal(q, N, true);
          });

        });
      } else {
        assert(!parent_locked); // TODO: need to prove that this assert can't trigger (i.e. need to prove bounding boxes become strictly smaller along any path)
        return insert_internal(q, T, false);
      }
    }
  }

  //assumes lock on T
  //uses q's interleave integer to recurse right or left
  pair<bool, bool> insert_internal(indexed_point q, node* T, bool parent_locked=false){
    node* N;
    if(lookup_bit(q.first, T->bit) == 0) N = T->Left();
    else N = T->Right();
    return insert_point0(q, T, N, T->bit-1, parent_locked);
  }

  bool insert_into_leaf(indexed_point q, node* parent, node* T){
      for(indexed_point p : T->Indexed_Pts()) {
        if(points_equal(p.second->pt, q.second->pt)){
          return false;
        }
      }

      std::vector<indexed_point> points;
      for(auto i : T->Indexed_Pts()) points.push_back(i);
      points.push_back(q);
      node* G = parent;
      bool left;
      if(G != nullptr) left = (G->Left() == T);
      //two cases: either (1) the leaf size is below the cutoff, in which case
      //we create a new leaf with the point q added, or (2) the leaf needs to be
      //split into an internal node and two leaf children
      if(T->size() + 1 < o_tree::node_cutoff || T->bit == 0){
        //case 1
        // std::cout << "regular insert" << std::endl;
        // std::cout << T->bit << std::endl;
        box b = T->Box();
        box bigger = box(b.first.minCoords(q.second->pt), b.second.maxCoords(q.second->pt));
        node* N = node::new_leaf(std::move(points), T->bit, bigger);
        assert(tree.load() == T || G != nullptr);
        if(G != nullptr) G->set_child(N, left);
        if(tree.load() == T) {
          // std::cout << "root changed (leaf: insert)" << std::endl; 
          set_root(N);
        }
        delete_single(T);
      } else {
        //sort points in leaf by interleave order
        int cut_point = 0;
        auto less_sort = [&] (indexed_point a, indexed_point b){
          return a.first < b.first;
        };
        std::sort(points.begin(), points.end(), less_sort);
        //search for cut point
        int new_bit = T->bit;
        while((cut_point == 0 | cut_point == points.size()) && new_bit != 0){
          size_t val = ((size_t) 1) << (new_bit - 1);
          size_t mask = (new_bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << new_bit);
          auto less = [&] (indexed_point x) {
            return (x.first & mask) < val;
          };
          cut_point = parlay::internal::binary_search(points, less);
          new_bit--;
        }
        std::vector<indexed_point> left_s;
        std::vector<indexed_point> right_s;
        for(size_t i=0; i<points.size(); i++){
          if(i<cut_point) left_s.push_back(points[i]);
          else right_s.push_back(points[i]);
        }
        node* L = node::new_leaf(std::move(left_s), new_bit);
        node* R = node::new_leaf(std::move(right_s), new_bit);
        node* P = node::new_node(L, R, new_bit+1);
        assert(tree.load() == T || G != nullptr);
        if(G != nullptr) G->set_child(P, left);
        if(T == tree.load()){
          set_root(P);
        }
        delete_single(T);
      }
      return true;
  }

  //should membership checking be part of deletion?
  //probably not bc can't differentiate between a deletion "in error"
  //vs a deletion that errors bc another concurrent process already deleted
  bool delete_point(vtx* p){
    return verlib::with_epoch([&] {
      int dims = p->pt.dimension();
      auto [b, Delta] = get_box_delta(dims);
      size_t interleave_int = o_tree::interleave_bits(p->pt, b.first, Delta);
      indexed_point q = std::make_pair(interleave_int, p);

      long cnt = 0;
      while(true) {
	if (cnt++ > 10000000) { std::cout << "overflow:" << parlay::worker_id() << ", " << tree.load()->is_leaf() << std::endl; cnt = 0;}
        node* R = tree.load();
        assert(R != nullptr);
        if(R->is_leaf()) {
          auto result = root_lock.try_lock_result([=] { 
              if(tree.load() != R) return -1; // validation failed
              else {
                bool success = R->lck.with_lock([=] { 
                  auto res = delete_from_leaf(q, nullptr, nullptr, R);
                  return res.first; 
                });
                if(success) return 1;
                else return 0;
              }
          });
          if(result.has_value() && result.value() != -1)
            return result.value() == 1;
        }
        else {
#ifdef PathCopy
          auto [a, locks_successful] = delete_point_path_copy(q, nullptr, R);
#else
          auto [a,b,locks_successful] = delete_point0(q, nullptr, R);
#endif
          if(locks_successful) return a;
        }
      }
    });
  }

  bool on_box_border(node* T, vtx* p){
    int dimensions = p->pt.dimension();
    auto box = T->Box();
    for (int i = 0; i < dimensions; i++) {
      if(box.first[i] == p->pt[i]) return true;
      if(box.second[i] == p->pt[i]) return true;
    }
    return false;
  }

  bool boxes_equal(box b1, box b2){
    int dims = b1.first.dimension();
    bool result = true;
    for(int i=0; i<dims; i++){
      result = result && (b1.first[i] == b2.first[i] && b1.second[i] == b2.second[i]);
    }
    return result;
  }

  //path copy version of delete
  pair<bool, bool> delete_point_path_copy(indexed_point q, node* parent, node* T){
    bool q_present; bool T_deleted = false; bool child_changed = true; bool locks_successful;
    node* N;

    if(lookup_bit(q.first, T->bit) == 0) N = T->Left();
    else N = T->Right();
    
    assert(!(T->is_leaf()));

    if(on_box_border(T, q.second) || N->is_leaf()) {
      // enter locking mode
      if(parent == nullptr) { // T is root
        // locks_taken_++;
        auto result = root_lock.try_lock_result([=] { 
          if(tree.load() != T) return make_pair(false, false); // validation failed
          else {
            std::vector<link> path; path.reserve(4);
            path.push_back((link){nullptr, false});
            bool success = delete_point_path_copy_helper(q, parent, T, path); 
            return make_pair(success, true);
          }
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); // locking failed 
      } else {
        // locks_taken_++;
        auto result = parent->lck.try_lock_result([=] {
          if(parent->removed.load() || (parent->Left() != T && parent->Right() != T)) { // if P is removed or P's child isn't T
            // std::cout << "failed valdiation" << std::endl; 
            return make_pair(false, false); // validation failed
          } else {
            std::vector<link> path; path.reserve(4);
            bool left_child = (parent->Left() == T);
            path.push_back((link){parent, left_child});
            bool success = delete_point_path_copy_helper(q, parent, T, path); 
            return make_pair(success, true);
          }
        });
        if(result.has_value()) return result.value();
        else return make_pair(false, false); // locking failed 
      }
    } else return delete_point_path_copy(q, T, N);
  }

  //returns pair of bools
  //first indicates whether q present in the data structure
  //second indicates whether the child node was changed--if not,
  //can end recursion early
  //third indicates if it failed to take all the locks
  bool delete_point_path_copy_helper(indexed_point q, node* parent, node* T, std::vector<link> &path){
    return T->lck.with_lock([=] {
      auto path_mutable = path;
      bool q_present; bool T_deleted = false; bool child_changed = true; bool locks_successful;
      node* N;
      bool go_left = lookup_bit(q.first, T->bit) == 0;
      if(go_left) N = T->Left();
      else N = T->Right();
      assert(!(T->is_leaf()));
      if(N->is_leaf()) // base case
        return delete_from_leaf_path_copy(q, T, parent, N, path_mutable);

      // recursive step
      link lnk = (link){T, go_left};
      path_mutable.push_back(lnk);
      return delete_point_path_copy_helper(q, T, N, path_mutable);
    });
  }

  //returns pair of bools
  //first indicates whether q present in the data structure
  //second indicates whether the child node was changed--if not,
  //can end recursion early
  //third indicates if it failed to take all the locks
  tuple<bool, bool, bool> delete_point0(indexed_point q, node* parent, node* T, bool parent_locked=false, bool T_locked=false){
    bool q_present; bool T_deleted = false; bool child_changed = true; bool locks_successful;
    node* N;

    if(lookup_bit(q.first, T->bit) == 0) N = T->Left();
    else N = T->Right();
    
    assert(!(T->is_leaf()));

    if(!parent_locked && (on_box_border(T, q.second) || N->is_leaf())) {
      // enter locking mode
      if(parent == nullptr) { // T is root
        // locks_taken_++;
        auto result = root_lock.try_lock_result([=] { 
          if(tree.load() != T) return make_tuple(false, false, false); // validation failed
          else return delete_point0(q, parent, T, true, false); 
        });
        if(result.has_value()) return result.value();
        else return make_tuple(false, false, false); // locking failed 
      } else {
        // locks_taken_++;
        auto result = parent->lck.try_lock_result([=] {
          if(parent->removed.load() || (parent->Left() != T && parent->Right() != T)) { // if P is removed or P's child isn't T
            // std::cout << "failed valdiation" << std::endl; 
            return make_tuple(false, false, false); // validation failed
          } else return delete_point0(q, parent, T, true, false); 
        });
        if(result.has_value()) return result.value();
        else return make_tuple(false, false, false); // locking failed 
      }
    }
    else if(parent_locked && !T_locked) {
      return T->lck.with_lock([=] {
        return delete_point0(q, parent, T, true, true); 
      });
    } else {
      if(N->is_leaf()) {
        auto [a,b] = delete_from_leaf(q, T, parent, N);
        q_present=a;
        T_deleted=b;
        locks_successful = true;
      } else{ 
        auto [a,b,c] =  delete_point0(q, T, N, parent_locked, false);
        q_present=a;
        child_changed=b;
        locks_successful = c;
        if(!locks_successful) return make_tuple(a,b,c);
      }

      //if the point wasn't present, no change needed
      //if child wasn't changed, no more changes needed, can end early
      if(!q_present || !child_changed) {return make_tuple(q_present, false, locks_successful);}
      
      if(!T_deleted) {
        //two cases:
        // (1) no change needed
        // (2) box needs to be made smaller
        if(on_box_border(T, q.second)){ //case 2
          // std::cout << "Internal delete, case 2" << std::endl;
          node* L = T->Left();
          node* R = T->Right();
          box b = box(L->Box().first.minCoords(R->Box().first),
            L->Box().second.maxCoords(R->Box().second));
          if(!boxes_equal(T->Box(), b)){
            node* N = node::new_node(L, R, T->bit);
            if(parent != nullptr){
              if(T == parent->Left()) parent->set_child(N, true);
              else parent->set_child(N, false);
            } else if(T == tree.load()){ 
              // std::cout << "root changed (internal, box change)" << std::endl;
              set_root(N);
            }
            delete_single(T);
          }
          return make_tuple(q_present, true, locks_successful);
        }
        else return make_tuple(q_present, false, locks_successful); //case 1
      } else return make_tuple(q_present, true, locks_successful);
    }
  }

  bool points_equal(point p, point q){
    int dims = p.dimension();
    bool result = true;
    for(int i=0; i<dims; i++){
      result = result && (p[i] == q[i]);
    }
    return result;
  }

  // grandparent is the last node in 'path'
  bool delete_from_leaf_path_copy(indexed_point q, node* parent, node* grandparent, node* T, std::vector<link> &path){
    std::vector<indexed_point> pts;
    bool cont = false;
    bool q_present = false;
    for(indexed_point p : T->Indexed_Pts()){
      if(cont) pts.push_back(p);
      else{
        if(points_equal(p.second->pt, q.second->pt)){
          q_present = true;
          cont = true;
        } else pts.push_back(p);    
      } 
    }
    if(!q_present) {return false;}
    node* sibling;
    int sib_size=0;
    bool sib_is_leaf=false;
    if(parent != nullptr){
      if(T == parent->Left()) sibling = parent->Right();
      else sibling = parent->Left();
      sib_size = sibling->size();
      if(sibling->is_leaf()) sib_is_leaf=true;
    }
    if(sib_is_leaf && sib_size + pts.size() < o_tree::node_cutoff_low && 
      parent != nullptr){
      for(indexed_point p : sibling->Indexed_Pts()){
        pts.push_back(p);
      }
      int bit;
      int dims = q.second->pt.dimension();
      if(grandparent != nullptr) bit = grandparent->bit-1;
      else bit = (o_tree::key_bits/dims)*dims;
      node* Leaf;
      if(on_box_border(parent, q.second)) Leaf = node::new_leaf(std::move(pts), bit);
      else Leaf = node::new_leaf(std::move(pts), bit, parent->Box());
      install_new_path(path, Leaf); // replace parent with Leaf
      delete_single(parent->Left()); // TODO: double check retires
      delete_single(parent->Right());
      delete_single(parent);
      return true;
    }
    else{
      if(pts.size() > 0){
        // std::cout << "Deletion from leaf, case 1" << std::endl;
        // node* Leaf = node::new_leaf(std::move(pts), T->bit);
        node* Leaf;
        if(on_box_border(T, q.second)) Leaf = node::new_leaf(std::move(pts), T->bit);
        else Leaf = node::new_leaf(std::move(pts), T->bit, T->Box());
        path.push_back((link){parent,parent->Left() == T});
        install_new_path(path, Leaf); // replace T with Leaf
        delete_single(T);
        return true;
      }
      else{
        // std::cout << "Deletion from leaf, case 3" << std::endl;
        if(parent == nullptr){
          std::cout << "ERROR: deleting entire tree not permitted" << std::endl;
          abort();
        }
        node* Leaf;
        if(sibling->is_leaf()){
          if(grandparent != nullptr){
            int bit = grandparent->bit-1;
            Leaf = node::new_leaf(sibling->Indexed_Pts(), bit, sibling->Box());
          }else{ //sibling becomes the root
            int dims = q.second->pt.dimension();
            int bit = dims*(o_tree::key_bits/dims);
            Leaf = node::new_leaf(sibling->Indexed_Pts(), bit, sibling->Box());
          }
          delete_single(sibling);
        } else Leaf = sibling;
        install_new_path(path, Leaf); // replace parent with Leaf
        // TODO: in some cases, sibling needs to be retired as well
        delete_single(parent);
        delete_single(T);
        return true;
      }
    } 
  }

  pair<bool, bool> delete_from_leaf(indexed_point q, node* parent, node* grandparent, node* T){
    std::vector<indexed_point> pts;
    bool cont = false;
    bool q_present = false;
    for(indexed_point p : T->Indexed_Pts()){
      if(cont) pts.push_back(p);
      else{
        if(points_equal(p.second->pt, q.second->pt)){
          q_present = true;
          cont = true;
        } else pts.push_back(p);    
      } 
    }
    if(!q_present) {return make_pair(false, false);}
    node* sibling;
    int sib_size=0;
    bool sib_is_leaf=false;
    if(parent != nullptr){
      if(T == parent->Left()) sibling = parent->Right();
      else sibling = parent->Left();
      sib_size = sibling->size();
      if(sibling->is_leaf()) sib_is_leaf=true;
    }
    if(sib_is_leaf && sib_size + pts.size() < o_tree::node_cutoff_low && 
      parent != nullptr){
      for(indexed_point p : sibling->Indexed_Pts()){
        pts.push_back(p);
      }
      int bit;
      int dims = q.second->pt.dimension();
      if(grandparent != nullptr) bit = grandparent->bit-1;
      else bit = (o_tree::key_bits/dims)*dims;
      node* Leaf;
      if(on_box_border(parent, q.second)) Leaf = node::new_leaf(std::move(pts), bit);
      else Leaf = node::new_leaf(std::move(pts), bit, parent->Box());
      if(grandparent != nullptr){
        if(parent == grandparent->Left()) grandparent->set_child(Leaf, true);
        else grandparent->set_child(Leaf, false);
      } else if(parent == tree.load()){ 
        // std::cout << "root changed (leaf)" << std::endl;
        set_root(Leaf);
      }
      delete_single(parent->Left());
      delete_single(parent->Right());
      delete_single(parent);
      return make_pair(q_present, true);
    }
    else{
      if(pts.size() > 0){
        // std::cout << "Deletion from leaf, case 1" << std::endl;
        // node* Leaf = node::new_leaf(std::move(pts), T->bit);
        node* Leaf;
        if(on_box_border(T, q.second)) Leaf = node::new_leaf(std::move(pts), T->bit);
        else Leaf = node::new_leaf(std::move(pts), T->bit, T->Box());
        if(parent != nullptr){
          if(T == parent->Left()) parent->set_child(Leaf, true);
          else parent->set_child(Leaf, false);
        } else if(T == tree.load()){ 
          // std::cout << "root changed (leaf was already root)" << std::endl;
          set_root(Leaf);
        }
        delete_single(T);
        return make_pair(q_present, false);
      }
      else{
        // std::cout << "Deletion from leaf, case 3" << std::endl;
        if(parent == nullptr){
          std::cout << "ERROR: deleting entire tree not permitted" << std::endl;
          abort();
        }
        if(sibling->is_leaf()){
          if(grandparent != nullptr){
            int bit = grandparent->bit-1;
            node* Leaf = node::new_leaf(sibling->Indexed_Pts(), bit, sibling->Box());
            if(parent == grandparent->Left()) grandparent->set_child(Leaf, true);
            else grandparent->set_child(Leaf, false);
          }else{ //sibling becomes the root
            int dims = q.second->pt.dimension();
            int bit = dims*(o_tree::key_bits/dims);
            node* Leaf = node::new_leaf(sibling->Indexed_Pts(), bit, sibling->Box());
            set_root(Leaf);
          }
          delete_single(sibling);
        } else{
          if(grandparent != nullptr){
            if(parent == grandparent->Left()) grandparent->set_child(sibling, true);
            else grandparent->set_child(sibling, false);
          } else{ 
            set_root(sibling);
          } 
        }
        delete_single(parent);
        delete_single(T);
        return make_pair(q_present, true);
      }
    }
  }

}; //this ends the k_nearest_neighbors structure

