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

  void set_box(box b){tree_box = b;}

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
    
    // looks for nearest neighbors for pt in Tree node T
    void k_nearest_rec(node* T) {
      box b;
      b = T->Box();
      // node::print_point(b.first);
      // node::print_point(b.second);
      // std::cout << max_distance << std::endl;
      if (within_epsilon_box(T, sqrt(max_distance))) { 
        // std::cout << "here" << std::endl;    
        if (report_stats) internal_cnt++;
        if (T->is_leaf()) {
          if (report_stats) leaf_cnt++;
          auto &Vtx = T->Indexed_Pts();
          for (int i = 0; i < T->size(); i++){
            // std::cout << T->size() << std::endl;
            // std::cout << Vtx.size() << std::endl;
            if (Vtx[i].second != vertex) { 
              // std::cout << "here" << std::endl;
              if (k < queue_cutoff){
                update_nearest(Vtx[i].second);
              } else{
                update_nearest_queue(Vtx[i].second);
            }
          } 
          }
          // abort();
        } else if (distance(T->Left()) < distance(T->Right())) {
          k_nearest_rec(T->Left());
          k_nearest_rec(T->Right());
        } else {
          k_nearest_rec(T->Right());
          k_nearest_rec(T->Left());
        }
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
      return std::move(range_candidates);
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
    box_delta bd = make_pair(b, Delta);
    return bd;
  }


  // takes in an integer and a position in said integer and returns whether the bit at that position is 0 or 1
  int lookup_bit(size_t interleave_integer, int pos){ //pos must be less than key_bits, can I throw error if not?
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
    // std::cout << Vtx[i].first << std::endl;
    std::cout << Vtx[i].second->identifier << std::endl;

  }
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
    verlib::with_snapshot([&] {
      kNN nn(p,k);
      nn.k_nearest_rec(tree.load()); //this is passing in a pointer to the o_tree root
      if (report_stats) {p->counter = nn.internal_cnt; p->counter2 = nn.leaf_cnt;}
      for (int i=0; i < k; i++)
        p->ngh[i] = nn[i];
    });
  }

  parlay::sequence<int> range_search(vtx* p, double rad){
      RNG rn(p,rad);
      rn.range_search_rec(tree.load()); //this is passing in a pointer to the o_tree root
      if (report_stats) {p->counter = rn.internal_cnt; p->counter2 = rn.leaf_cnt;}
      return std::move(rn.return_answer());
  }
  

 
  parlay::sequence<vtx*> z_sort(parlay::sequence<vtx*> v, box b){ 
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
    return std::move(v3); 
  }

  using indexed_point = typename o_tree::indexed_point; 

  bool within_box(node* T, vtx* vertex) {
      int dimensions = vertex->pt.dimension();
      auto box = T->Box();
      bool result = true;
      for (int i = 0; i < dimensions; i++) {
        result = (result &&
          (box.first[i] <= vertex->pt[i]) &&
          (box.second[i] >= vertex->pt[i]));
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

  //takes in a single point and a pointer to the root of the tree
  //as well as a bounding box and its largest side length
  //assumes integer coordinates
  void insert_point(vtx* p){
    verlib::with_epoch([&] {
      int dims = p->pt.dimension();
      auto [b, Delta] = get_box_delta(dims);
      size_t interleave_int = o_tree::interleave_bits(p->pt, b.first, Delta);
      indexed_point q = std::make_pair(interleave_int, p);
      while(true) {
        node* R = tree.load();
        if(insert_point0(q, nullptr, R, o_tree::key_bits, false)) break;
        // std::cout << "insert attempt failed" << std::endl; 
      }
      // while(!insert_point0(q, nullptr, R, o_tree::key_bits, false)) {} // repeat until success
    });
  }

  bool insert_point0(indexed_point q, node* parent, node* T, int bit, bool parent_locked=false){
    //lock T
    // assert(false);
    if(!parent_locked && (T->is_leaf() || !within_box(T, q.second))) {
      // enter locking mode
      if(parent == nullptr) { // T is root
        return root_lock.try_lock([=] { 
          if(tree.load() != T) return false;
          else return insert_point0(q, parent, T, bit, true); 
        });
      } else {
        return parent->lck.try_lock([=] {
          if(parent->removed.load() || !(parent->Left() == T || parent->Right() == T)) { // if P is removed or P's child isn't T
            // std::cout << "failed valdiation" << std::endl; 
            return false;
          } else return insert_point0(q, parent, T, bit, true);
        });
      }
    }

    if(T->is_leaf()) {
      assert(parent_locked);
      return T->lck.with_lock([=] { return insert_into_leaf(q, parent, T); });
    } else {
      if(!within_box(T, q.second)) {
        assert(parent_locked);
        return T->lck.with_lock([=] {
          //two cases: (1) need to create new leaf and internal node,
          //or (2) bit unchanged and box changed

          //if next bit of integer isn't same as node,
          //iterate until either make a leaf or bits match
          node* Q = T;
          while(!(Q->is_leaf())){node* L = Q->Right(); Q=L;}
          indexed_point s = Q->Indexed_Pts()[0];
          int cur_bit = bit;
          assert(cur_bit >= T->bit);
          if(parent != nullptr) assert(cur_bit < G->bit);
          while(cur_bit > T->bit) {
            if(lookup_bit(q.first, cur_bit) != lookup_bit(s.first, cur_bit)){
              // std::cout << "case 1" << std::endl;
              //we know we are in case 1
              //form leaf
              node* G = parent;
              parlay::sequence<indexed_point> points = {q};
              node* R = node::new_leaf(parlay::make_slice(points), cur_bit-1);
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
              return true;
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
      } else return insert_internal(q, T, false);
    }
  }

  //assumes lock on T
  //uses q's interleave integer to recurse right or left
  bool insert_internal(indexed_point q, node* T, bool parent_locked=false){
    node* N;
    if(lookup_bit(q.first, T->bit) == 0) N = T->Left();
    else N = T->Right();
    return insert_point0(q, T, N, T->bit-1, parent_locked);
  }

  bool insert_into_leaf(indexed_point q, node* parent, node* T){
      auto points = parlay::tabulate(T->size()+1, [&] (long i) {
	  return (i == T->size()) ? q : T->Indexed_Pts()[i];}, 1000);
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
        node* N = node::new_leaf(std::move(points), T->bit, std::move(bigger));
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
        parlay::slice<indexed_point*, indexed_point*> pts = parlay::make_slice(points);
        auto left_s = parlay::tabulate(cut_point, [&] (size_t i) {return points[i];}, 1000);
        auto right_s = parlay::tabulate(points.size()-cut_point, [&] (size_t i) {return points[cut_point+i];}, 1000);
        node* L = node::new_leaf(std::move(left_s), new_bit);
        node* R = node::new_leaf(std::move(right_s), new_bit);
        node* P = node::new_node(L, R, new_bit+1);
        assert(tree.load() == T || G != nullptr);
        if(G != nullptr) G->set_child(P, left);
        if(T == tree.load()){
          // std::cout << "root changed (leaf: split and insert)" << std::endl; 
          set_root(P);
        }
        delete_single(T);
      }
      return true;
  }

  //should membership checking be part of deletion?
  //probably not bc can't differentiate between a deletion "in error"
  //vs a deletion that errors bc another concurrent process already deleted
  void delete_point(vtx* p){
    int dims = p->pt.dimension();
    auto [b, Delta] = get_box_delta(dims);
    size_t interleave_int = o_tree::interleave_bits(p->pt, b.first, Delta);
    indexed_point q = std::make_pair(interleave_int, p);
    node* R = tree.load();
    assert(R != nullptr);
    if(R->is_leaf()) {
      delete_from_leaf(q, nullptr, nullptr, R);
    }
    else delete_point0(q, nullptr, nullptr, R);
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

  //returns pair of bools
  //first indicates whether q present in the data structure
  //second indicates whether the child node was changed--if not,
  //can end recursion early
  std::pair<bool, bool> delete_point0(indexed_point q, node* parent, node* grandparent, node* T){
    bool q_present; bool T_deleted = false; bool child_changed = true;
    node* N;
    if(lookup_bit(q.first, T->bit) == 0) N = T->Left();
    else N = T->Right();
    if(N->is_leaf()) {
      auto [a,b] = delete_from_leaf(q, T, parent, N);
      q_present=a;
      T_deleted=b;
    } else{ 
      auto [a,b] =  delete_point0(q, T, parent, N);
      q_present=a;
      child_changed=b;
    }

    //if the point wasn't present, no change needed
    //if child wasn't changed, no more changes needed, can end early
    if(!q_present || !child_changed) {return std::make_pair(q_present, false);}
    
    if(!T_deleted){
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
        return std::make_pair(q_present, true);
      }
      else return std::make_pair(q_present, false); //case 1
    } else return std::make_pair(q_present, true);
  }

  bool points_equal(point p, point q){
    int dims = p.dimension();
    bool result = true;
    for(int i=0; i<dims; i++){
      result = result && (p[i] == q[i]);
    }
    return result;
  }

  std::pair<bool, bool> delete_from_leaf(indexed_point q, node* parent, node* grandparent, node* T){
    parlay::sequence<indexed_point> pts;
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
    if(!q_present) {return std::make_pair(false, false);}
    node* sibling;
    int sib_size=0;
    bool sib_is_leaf=false;
    if(parent != nullptr){
      if(T == parent->Left()) sibling = parent->Right();
      else sibling = parent->Left();
      sib_size = sibling->size();
      if(sibling->is_leaf()) sib_is_leaf=true;
    }
    if(sib_is_leaf && sib_size + pts.size() < o_tree::node_cutoff && 
      parent != nullptr){
      // std::cout << "Deletion from leaf, case 2" << std::endl;
      for(indexed_point p : sibling->Indexed_Pts()){
        pts.push_back(p);
      }
      int bit;
      int dims = q.second->pt.dimension();
      if(grandparent != nullptr) bit = grandparent->bit-1;
      else bit = (o_tree::key_bits/dims)*dims;
      node* Leaf;
      if(on_box_border(parent, q.second)) Leaf = node::new_leaf(std::move(pts), bit);
      else Leaf = node::new_leaf(std::move(pts), bit, std::move(parent->Box()));
      // node* Leaf = node::new_leaf(std::move(pts), bit);
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
      return std::make_pair(q_present, true);
    }
    else{
      if(pts.size() > 0){
        // std::cout << "Deletion from leaf, case 1" << std::endl;
        // node* Leaf = node::new_leaf(std::move(pts), T->bit);
        node* Leaf;
        if(on_box_border(T, q.second)) Leaf = node::new_leaf(std::move(pts), T->bit);
        else Leaf = node::new_leaf(std::move(pts), T->bit, std::move(T->Box()));
        if(parent != nullptr){
          if(T == parent->Left()) parent->set_child(Leaf, true);
          else parent->set_child(Leaf, false);
        } else if(T == tree.load()){ 
          // std::cout << "root changed (leaf was already root)" << std::endl;
          set_root(Leaf);
        }
        delete_single(T);
        return std::make_pair(q_present, false);
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
            // node* Leaf = node::new_leaf(parlay::make_slice(sibling->Indexed_Pts()), bit);
            // node* Leaf = node::new_leaf(parlay::make_slice(sibling->Indexed_Pts()), bit);
            // node* Leaf = node::new_leaf(parlay::make_slice(sibling->Indexed_Pts()), bit, sibling->Box());
            node* Leaf = node::new_leaf(std::move(sibling->Indexed_Pts()), bit, std::move(sibling->Box()));
            if(parent == grandparent->Left()) grandparent->set_child(Leaf, true);
            else grandparent->set_child(Leaf, false);
          }else{ //sibling becomes the root
            int dims = q.second->pt.dimension();
            int bit = dims*(o_tree::key_bits/dims);
            node* Leaf = node::new_leaf(std::move(sibling->Indexed_Pts()), bit, std::move(sibling->Box()));
            // std::cout << "Root changed: sibling promoted to root" << std::endl;
            set_root(Leaf);
          }
          delete_single(sibling);
        } else{
          if(grandparent != nullptr){
            if(parent == grandparent->Left()) grandparent->set_child(sibling, true);
            else grandparent->set_child(sibling, false);
          } else{ 
            // std::cout << "Root changed: sibling promoted to root" << std::endl;
            set_root(sibling);
          } 
        }
        delete_single(parent);
        delete_single(T);
        return std::make_pair(q_present, true);
      }
    }
  }

}; //this ends the k_nearest_neighbors structure

