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
#include <math.h>  
#include <mutex>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/get_time.h"
#include "common/atomics.h"
#include "parlay/hash_table.h"

// vtx must support v->pt
// and v->pt must support pt.dimension(), pt[i],
//    (pt1 - pt2).Length(), pt1 + (pt2 - pt3)
//    pt1.minCoords(pt2), pt1.maxCoords(pt2),


template <typename vtx>  
struct ck_tree{

	using point = typename vtx::pointT;
  	using box = std::pair<point,point>;
  	using slice_t = decltype(make_slice(parlay::sequence<vtx*>()));


	static box get_box(slice_t V) { // parlay::sequence<vtx*> &V) {
		size_t n = V.size();
		auto minmax = [&] (box x, box y) {
		  return box(x.first.minCoords(y.first),
		 x.second.maxCoords(y.second));};

		// uses a delayed sequence to avoid making a copy
		auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
		return box(V[i]->pt, V[i]->pt);});
		box identity = pts[0];
		return parlay::reduce(pts, parlay::make_monoid(minmax,identity));
	}

	struct node{


		// struct hashNode{
		// 	using eType = node*;
		// 	using kType = node*; 
		// 	eType empty() {return NULL;}
		// 	kType getKey(eType v) {return v;}
		// 	size_t hash(kType s) {return parlay::hash64((s->center())[0]);} //+(s->second->center())[0]);}
		// 	bool replaceQ(eType, eType) { return 0; }
	 //  		bool cas(eType* p, eType o, eType n) {
		// 	    return pbbs::atomic_compare_and_swap(p, o, n);
		// 	}
		// 	int cmp(kType s, kType s2) {
		// 	    return ((s->center())[0] > (s2->center())[0]) ? 1 : (((s->center())[0] == (s2->center())[0]) ? 0 : -1);
		// 	}
		// };

		// typedef parlay::hashtable<hashNode> NodeTable;
		// static NodeTable makeNodeTable(size_t m) {
		//  	return NodeTable(m, hashNode());
		// }

		public:
		    box Box() {return b;}
		    point center() {return centerv;}
		    double Max_dim() {return max_dim;}
		    double Diameter() {return diameter;}
		    size_t size() {return n;}
		    size_t ID() {return id;}
		    bool is_leaf() {return L == nullptr;}
		    node* Left() {return L;}
		    node* Right() {return R;}
		    node* Parent() {return parent;}
		    vtx* Vertex() {return V;}
			// interactions_wrapper I; 
			parlay::sequence<node*> Interactions() {return interactions;}
			// std::mutex m; 
			// slice_n Interactions() {return node_table.entries();}

		void add_interaction(node* T){
			// std::mutex m; 
			// std::scoped_lock lock(m);
			interactions.push_back(T);
		}

			// void add_interaction(node* T){
			// 	node_table.insert(T);
			// }


		// construct a leaf node with a sequence of points directly in it
	    node(slice_t Pts, size_t idty) { 
	    	n = 1;
	      	parent = nullptr;
	     	L = R = nullptr;
	      	b = get_box(Pts);
	      	set_center();
	      	set_max_dim();
	      	set_diameter();
	      	interactions = parlay::sequence<node*>();
	      	// I = interactions_wrapper();
	      	// node_table = makeNodeTable(100);
	      	V = Pts[0];
	      	id = idty;
	    }

	    // construct an internal binary node
	    node(node* L, node* R, size_t idty) : L(L), R(R) { 
		    parent = nullptr;
		    b = box(L->b.first.minCoords(R->b.first),
			  	L->b.second.maxCoords(R->b.second));
		    n = L->size() + R->size();
		    set_center();
		    set_max_dim(); 
		    set_diameter();
		    interactions = parlay::sequence<node*>();
		    // I = interactions_wrapper(); 
		    // node_table = makeNodeTable(100);
		    id = idty; 
	    }


	    static node* new_leaf(slice_t Pts, size_t ID) {
      		node* r = alloc_node();
      		new (r) node(Pts, ID);
      		return r;
    	}

    	static node* new_node(node* L, node* R, size_t ID) {
      		node* nd = alloc_node();
      		new (nd) node(L, R, ID);
      		// both children point to this node as their parent
      		L->parent = R->parent = nd;
     		return nd;
    	}

    	~node() {
      		// need to collect in parallel
      		parlay::par_do_if(n > 1000,
			[&] () { delete_tree(L);},
			[&] () { delete_tree(R);});
    	}

    	// map a function f(p,node_ptr) over the points, passing
	    // in a pointer to a vertex, and a pointer to the leaf node it is in.
	    // f should return void

	    //pass in a function to compute nearest neighbors
	    template <typename F>
	    void map(F f) { 
	      	if (is_leaf()){ 
	      		// std::cout << "at map" << std::endl; 
	      		f(V, this);
	      	}
	      	else {
				parlay::par_do_if(n > 10,
					[&] () {L->map(f);},
					[&] () {R->map(f);});
	      	}
	    }

    	size_t depth() {
      		if (is_leaf()) return 0;
      		else {
				size_t l, r;
				parlay::par_do_if(n > 1000,
			  		[&] () {l = L->depth();},
			  		[&] () {r = R->depth();});
			return 1 + std::max(l,r);
      		}
    	}

    	// recursively frees the tree
    	static void delete_tree(node* T) {
      		if (T != nullptr) {
				T->~node();
				node::free_node(T);
      		}
    	}		

	    // disable copy and move constructors/assignment since
	    // they are dangerous with with free.
	    node(const node&) = delete;
	    node(node&&) = delete;
	    node& operator=(node const&) = delete;
	    node& operator=(node&&) = delete;


	private:

    	size_t n;
    	double diameter; 
    	double max_dim; 
    	size_t id; 
    	node *parent;
    	node *L;
    	node *R;
    	box b;
    	point centerv;
    	vtx* V; 
    	parlay::sequence<node*> interactions;
    	// NodeTable node_table; 

    	void set_center() {			   
      		centerv = b.first + (b.second-b.first)/2;
    	}

    	void set_max_dim(){
			double maxdim = 0; 
			int dim = b.first.dimension();
			for(int i=0; i<dim; i++){
				double new_max = abs(b.first[i]-b.second[i]);
				if (new_max > maxdim){
					maxdim = new_max;
				}
			}
			max_dim = maxdim; 
		}

		void set_diameter(){
			double diam = (b.first-b.second).Length();
			diameter = diam; 
		}

   		static parlay::type_allocator<node> node_allocator;
    	static node* alloc_node() {
      		return node_allocator.alloc();}

    	static void free_node(node* T) {
      		node_allocator.free(T);}

	}; //ends node structure

	static bool well_separated(node* A, node* B, double s){
		double m_a = A->Diameter();
		double m_b = B->Diameter();
		double diameter = max(m_a, m_b); //diameter of the smallest sphere that can capture each box
		double d = (A->center()- B->center()).Length(); //distance between the centers of the two boxes
		return (d-diameter >= .5*s*diameter); //is the distance between the two balls larger than .5*s*r
	}

	// A unique pointer to a tree node to ensure the tree is
  	// destructed when the pointer is, and that no copies are made.
  	struct delete_tree {void operator() (node *T) const {node::delete_tree(T);}};
  	using tree_ptr = std::unique_ptr<node,delete_tree>;

 //  	using slice_n = decltype(make_slice(parlay::sequence<node*>()));
 //  	using node_pair = std::pair<node*, node*>;

	// struct hashNode{
	// 	using eType = node_pair*;
	// 	using kType = node_pair*; 
	// 	eType empty() {return NULL;}
	// 	kType getKey(eType v) {return v;}
	// 	size_t hash(kType s) {return parlay::hash64(rand());}
	// 	// size_t hash(kType s) {return parlay::hash64((s->first->center())[0]+(s->second->center())[0]+rand());}
	// 	// size_t hash(kType s) {return std::hash<double>{}(31*std::hash<double>{}(31*std::hash<double>{}(23)
	// 	// 	+(s->first->center())[0])+(s->second->center())[1]);}
	// 	bool replaceQ(eType, eType) { return 0; }
 //  		bool cas(eType* p, eType o, eType n) {
	// 	    return pbbs::atomic_compare_and_swap(p, o, n);
	// 	}
	// 	int cmp(kType s, kType s2) {
	// 	    // return ((s->first->ID()+3*(s->second->ID()) > s2->first->ID()+3*(s2->second->ID())) ? 1
	// 	    //  : (s->first->ID()+3*(s->second->ID())== s2->first->ID()+3*(s2->second->ID())) ? 0 : -1);
	// 	    return 1; 
	// 	}
	// };

	// typedef parlay::hashtable<hashNode> NodeTable;
	// static NodeTable makeNodeTable(size_t m) {
	//  	return NodeTable(m, hashNode());
	// }

	// static void ints_rec(parlay::slice<node_pair*> ints){
	// 	return; 
	// }

	// static void apply_interactions(node* T, NodeTable interactionsTable){
	// 	std::cout << interactionsTable.count() << std::endl; 
	// 	auto less = [&] (node_pair* a, node_pair* b){
	// 		return a->first->ID() < b->first->ID();
	// 	};
	// 	// std::cout << "made lambda" << std::endl; 
	// 	// auto x = parlay::sort(interactionsTable.entries(), less);
	// 	// std::cout << "here" << std::endl; 
	// 	// std::cout << x[0]->first->Vertex()->pt[0];
	// 	// ints_rec(parlay::make_slice(x));
	// }

	// static void wsr_wrapper(size_t n, node* T, double s, int k){
	// 	NodeTable interactionsTable = makeNodeTable(10*k*n);
	// 	wsr(T, s, k, interactionsTable);
	// 	std::cout << "done with interactions" << std::endl; 
	// 	apply_interactions(T, interactionsTable);
	// }


  	static void wsr(node* T, double s, int k){ //, NodeTable &interactionsTable){
  		// NodeTable table = makeNodeTable(1000);
  		if(T->is_leaf()){
  			return;
  		} else{
			node* L = T -> Left();
			node* R = T -> Right();
			wsrChildren(L, R, s, k);
			size_t n = T->size();
			parlay::par_do_if(n > 100, 
				[&] () {wsr(L, s, k);},
				[&] () {wsr(R, s, k);}
			);
  		}
	}

	//s should be set greater than 2 for nearest neighbor applications
	static void wsrChildren(node* L, node* R, double s, int k){ //, NodeTable &interactionsTable){ 	
		if(well_separated(L, R, s)){
			if(L->size() <= k){
				// std::scoped_lock lock(L->m);
				// std::cout << "count before first insertion " << interactionsTable.count() << std::endl; 
				// node_pair pair = std::make_pair(L, R);
				// node_pair* pair_ptr = &pair;
				// interactionsTable.insert(pair_ptr);
				// std::cout << "count after first insertion " << interactionsTable.count() << std::endl; 
				L -> add_interaction(R);
			}
			if(R->size() <= k){
				// std::scoped_lock lock(R->m);
				// node_pair pair = std::make_pair(R, L);
				// node_pair* pair_ptr = &pair;
				// interactionsTable.insert(pair_ptr);
				R -> add_interaction(L);
			}
		} else{
			size_t n = L->size() + R->size();
			double m_L = L->Max_dim();
			double m_R = R->Max_dim();
			if(m_L > m_R){
				// parlay::par_do_if(n > 100, 
				// 	[&] () {wsrChildren(R, L->Left(), s);},
				// 	[&] () {wsrChildren(R, L->Right(), s);}
				// );
	  			wsrChildren(R, L->Left(), s, k);
	  			wsrChildren(R, L->Right(), s, k);
			} else{
				// parlay::par_do_if(n > 100, 
				// 	[&] () {wsrChildren(L, R->Left(), s);},
				// 	[&] () {wsrChildren(L, R->Right(), s);}
				// );
	  			wsrChildren(L, R->Left(), s, k);
	  			wsrChildren(L, R->Right(), s, k);
			}
		}
	}

	// build a tree given a sequence of pointers to points
  	template <typename Seq>
  	static tree_ptr build(Seq &P) {
    	timer t("oct_tree", false);
    	size_t n = P.size();
  		parlay::sequence<vtx*> Tmp;
  		Tmp = parlay::sequence<vtx*>(n);
    	node* r = build_recursive(parlay::make_slice(P), parlay::make_slice(Tmp), 0);
    	t.next("build");
    	return tree_ptr(r);
  	}

	static bool are_equal(point p, point q){
		int d = p.dimension();
		for (int j = 0; j<d; j++){
			if(p[j] != q[j]){
				return false;
			}
		}
		return true;
	}

  	static node* build_recursive(slice_t points, slice_t Tmp, size_t id_offset){
  		if(points.size() == 0){
  			std::cout << "ERROR: passed in slice of size 0 when building tree" << std::endl; 
  			abort(); 
  		}
		if(points.size() == 1){
  			return node::new_leaf(points, id_offset);
  		}else{
	  		//identify the box around the points and its largest dimension
	  		int dim = (points[0]->pt).dimension();
	  		box b = get_box(points);
	  		size_t d = 0;
	  		double Delta = 0.0;
	  		for (int i=0; i < dim; i++) {
	    		if (std::abs(b.second[i] - b.first[i]) > Delta) {
	      			d = i;
	      			Delta = std::abs(b.second[i] - b.first[i]);
	    		}
	  		}
	  		double splitpoint = (b.first[d] + b.second[d])/2.0;
	  		//pack points into new arrays based on splitpoint
	  		auto flagsLeft = parlay::map(points, [&] (vtx* p) -> bool {
	      		return p->pt[d] < splitpoint;});
	  		auto flagsRight = parlay::delayed_map(flagsLeft, [] (bool x) {
	      		return !x;});
	  		size_t n = points.size();
	  		size_t nl = parlay::pack_into(points, flagsLeft, Tmp);
	  		parlay::pack_into(points, flagsRight, Tmp.cut(nl, n));
	  		if(nl == 0){ //edge case where all points are the same
	  			nl = n/2; 
	  		}
	  		//create left and right children
  			node *L, *R; 
	  		parlay::par_do_if(n > 1000, 
	  			[&] () {L = build_recursive(Tmp.cut(0, nl), points.cut(0, nl), id_offset);},
	  			[&] () {R = build_recursive(Tmp.cut(nl, n), points.cut(nl, n), id_offset+nl);}
	  		);
	  		//create parent node
	  		return node::new_node(L, R, (nl+id_offset)*2);
	  	}
  	}

}; //ends CKtree structure





  	// 	// build a tree given a sequence of pointers to points
  	// template <typename Seq>
  	// static tree_ptr build_ck_0(Seq &P) {
   //  	timer t("oct_tree", false);
   //  	size_t n = P.size();
  	// 	parlay::sequence<vtx*> Tmp;
  	// 	Tmp = parlay::sequence<vtx*>(n);
   //  	node* r = build_recursive_ck_0(parlay::make_slice(P), parlay::make_slice(Tmp));
   //  	t.next("build");
   //  	return tree_ptr(r);
  	// }

  	// static node* build_recursive_ck_0(slice_t points, slice_t Tmp){
  	// 	size_t initial_size = points.size();
  	// 	size_t current_size = initial_size; 
  	// 	while(current_size > initial_size/2){
  	// 		int dim = (points[0]->pt).dimension();
	  // 		box b = get_box(points);
	  // 		size_t d = 0;
	  // 		double Delta = 0.0;
	  // 		for (int i=0; i < dim; i++) {
	  //   		if (std::abs(b.second[i] - b.first[i]) > Delta) {
	  //     			d = i;
	  //     			Delta = std::abs(b.second[i] - b.first[i]);
	  //   		}
	  // 		}
	  // 		// std::cout << d << std::endl; 
	  // 		double splitpoint = (b.first[d] + b.second[d])/2.0;
	  // 		// std::cout << splitpoint << std::endl; 
	  // 		//pack points into new arrays based on splitpoint
	  // 		auto flagsLeft = parlay::map(points, [&] (vtx* p) -> bool {
	  //     		return p->pt[d] < splitpoint;});
	  // 		auto flagsRight = parlay::delayed_map(flagsLeft, [] (bool x) {
	  //     		return !x;});
	  // 		size_t n = points.size();
	  // 		size_t nl = parlay::pack_into(points, flagsLeft, Tmp);
	  // 		parlay::pack_into(points, flagsRight, Tmp.cut(nl, n));
	  // 		if(nl == 0){ //edge case where all points are the same
	  // 			nl = n/2; 
	  // 		}
  	// 	}
  	// }

  	// // build a tree given a sequence of pointers to points
  	// template <typename Seq>
  	// static tree_ptr build_ck(Seq &P) {
   //  	timer t("ck_tree", false);
   //  	size_t n = P.size();
   //  	int dim = P[0]->pt.dimension();
   //  	parlay::sequence<slice_t> slice_seq;
   //  	slice_seq = parlay::sequence<slice_t>(dim);
   //  	parlay::sequence<slice_t> slice_seq_tmp;
   //  	slice_seq_tmp = parlay::sequence<slice_t>(dim);
   //  	for(int i=0; i<dim; i++){
   //  		slice_seq[i] = sort_coord(P, i);
   //  		parlay::sequence<vtx*> Tmp;
  	// 		Tmp = parlay::sequence<vtx*>(n);
  	// 		slice_seq_tmp[i] = Tmp; 
   //  	}
   //  	node* r = build_ck_rec(slice_seq, slice_seq_tmp, dim);
   //  	t.next("build");
   //  	return tree_ptr(r);
  	// }

  	// slice_t sort_coord(parlay::sequence<vtx*> v, int dim){
  	// 	auto less = [] (vtx* a, vtx* b){
  	// 		return (a->pt)[dim] < (b->pt)[dim];
  	// 	};
  	// 	auto x = parlay::sort(v, less);
  	// 	return parlay::make_slice(x);
  	// }

  	// static node* build_ck_rec(parlay::sequence<slice_t> slice_seq, parlay::sequence<slice_t> size_seq_tmp, int dims){
  	// 	size_t size_init = slice_seq[0].size();
  	// 	size_t running_size = size_init; 
  	// 	while(running_size > size_init/2){
  	// 		double maxdim = 0; 
  	// 		double d = 0; 
  	// 		for(int i=0; i<dims; i++){

  	// 		}
  	// 		double splitpoint = maxdim/2;
  	// 		auto less = [&] (vtx* a){
  	// 			return (a->pt)[d] < splitpoint; 
  	// 		};
  	// 		size_t splitcoord = parlay::internal::binary_search(slice_seq[d], less);
  	// 		//create new wrappers for the sorted sequences
  	// 		parlay::sequence<slice_t> slice_seq_less = parlay::sequence<slice_t>(dims);
  	// 		parlay::sequence<slice_t> slice_seq_tmp_less = parlay::sequence<slice_t>(dims);
  	// 		parlay::sequence<slice_t> slice_seq_greater = parlay::sequence<slice_t>(dims);
  	// 		parlay::sequence<slice_t> slice_seq_tmp_greater = parlay::sequence<slice_t>(dims);
  	// 		//split the sorted sequence we cut on
  	// 		slice_seq_less[d] = slice_seq[d].cut(0, splitcoord);
  	// 		slice_seq_tmp_less[d] = slice_seq_tmp[d].cut(0, splitcoord);
  	// 		slice_seq_greater[d] = slice_seq[d].cut(splitcoord, running_size);
  	// 		slice_seq_tmp_less[d] = slice_seq_tmp[d].cut(splitcoord, running_size);
  	// 		for(int i=0; i<dims; i++){
  	// 			if(i != d){

  	// 			}
  	// 		}
  	// 	}
  	// }


		// parlay::hashtable<hash_node> table; 

		// // Example for hashing numeric values.
		// // T must be some integer type
		// template <class T>
		// struct hash_numeric {
		//   using eType = T;
		//   using kType = T;
		//   eType empty() { return -1; }
		//   kType getKey(eType v) { return v; }
		//   size_t hash(kType v) { return static_cast<size_t>(hash64(v)); }
		//   int cmp(kType v, kType b) { return (v > b) ? 1 : ((v == b) ? 0 : -1); }
		//   bool replaceQ(eType, eType) { return 0; }
		//   eType update(eType v, eType) { return v; }
		//   bool cas(eType* p, eType o, eType n) {
		//     // TODO: Make this use atomics properly. This is a quick
		//     // fix to get around the fact that the hashtable does
		//     // not use atomics. This will break for types that
		//     // do not inline perfectly inside a std::atomic (i.e.,
		//     // any type that the standard library chooses to lock)
		//     return std::atomic_compare_exchange_strong_explicit(
		//       reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
		//   }
		// };

		// parlay::hashtable<hash_numeric<int>> table;

		// struct interactions_wrapper{
		// 	parlay::sequence<node*> interactions;
		// 	std::mutex *m; 

		// 	interactions_wrapper(){
		// 		interactions = parlay::sequence<node*>();
		// 		m = new std::mutex();
		// 	}


		// 	void add_interaction(node* T){
		// 		std::scoped_lock lock(*m);
		// 		interactions.push_back(T);
		// 	}

		// };


