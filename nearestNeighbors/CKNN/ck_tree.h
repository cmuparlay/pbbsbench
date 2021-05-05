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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/get_time.h"

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

		public:
		    box Box() {return b;}
		    point center() {return centerv;}
		    double Max_dim() {return max_dim;}
		    double Diameter() {return diameter;}
		    size_t size() {return n;}
		    bool is_leaf() {return L == nullptr;}
		    node* Left() {return L;}
		    node* Right() {return R;}
		    node* Parent() {return parent;}
		    parlay::sequence<node*> Interactions() {return interactions;}
		    vtx* Vertex() {return V;}

		void add_interaction(node* T){
			size_t len = interactions.size();
			for(size_t i=0; i<len; i++){
				if(interactions[i] == T){
					return;
				}
			}
			parlay::sequence<node*> ints;
			ints = parlay::sequence<node*>(len+1);
			for(size_t i=0; i<len; i++){
				ints[i] = interactions[i];
			}
			ints[len] = T;
			interactions = ints; 
		}

		// construct a leaf node with a sequence of points directly in it
	    node(slice_t Pts) { 
	    	n = 1;
	      	parent = nullptr;
	     	L = R = nullptr;
	      	b = get_box(Pts);
	      	set_center();
	      	set_max_dim();
	      	set_diameter();
	      	interactions = parlay::sequence<node*>();
	      	V = Pts[0];
	    }

	    // construct an internal binary node
	    node(node* L, node* R) : L(L), R(R) { 
		    parent = nullptr;
		    b = box(L->b.first.minCoords(R->b.first),
			  	L->b.second.maxCoords(R->b.second));
		    n = L->size() + R->size();
		    set_center();
		    set_max_dim(); 
		    set_diameter();
		    interactions = parlay::sequence<node*>();
	    }

	    static node* new_leaf(slice_t Pts) {
      		node* r = alloc_node();
      		new (r) node(Pts);
      		return r;
    	}

    	static node* new_node(node* L, node* R) {
      		node* nd = alloc_node();
      		new (nd) node(L, R);
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
	      	if (is_leaf()) f(V, this);
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
    	node *parent;
    	node *L;
    	node *R;
    	box b;
    	point centerv;
    	vtx* V; 
    	parlay::sequence<node*> interactions;

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

    	// uses the parlay memory manager, could be replaced will alloc/free
   		static parlay::type_allocator<node> node_allocator;
    	static node* alloc_node() {
      		return node_allocator.alloc();}

    	static void free_node(node* T) {
      		node_allocator.free(T);}

	}; //ends node structure

	static bool well_separated(node* A, node* B, double s){
		// std::cout << "here1" << std::endl; 
		double m_a = A->Diameter();
		// std::cout << "here2" << std::endl; 
		// std::cout << B->is_leaf() << std::endl; 
		double m_b = B->Diameter();
		// std::cout << "here3" << std::endl; 
		double diameter = max(m_a, m_b); //diameter of the smallest sphere that can capture each box
		// std::cout << "here4" << std::endl; 
		double d = (A->center()- B->center()).Length(); //distance between the centers of the two boxes
		// std::cout << "here5" << std::endl; 
		return (d-diameter >= .5*s*diameter); //is the distance between the two balls larger than .5*s*r
	}

	// A unique pointer to a tree node to ensure the tree is
  	// destructed when the pointer is, and that no copies are made.
  	struct delete_tree {void operator() (node *T) const {node::delete_tree(T);}};
  	using tree_ptr = std::unique_ptr<node,delete_tree>;

  	static void wsr(node* T, double s){
  		if(T->is_leaf()){
  			return;
  		} else{
  			size_t n = T->size();
			node* L = T -> Left();
			node* R = T -> Right();
			wsrChildren(L, R, s);
			parlay::par_do_if(n > 1000, 
				[&] () {wsr(L, s);},
				[&] () {wsr(R, s);}
			);
  		}
	}

	static void wsrChildren(node* L, node* R, double s){ //s should be set greater than 2 for nearest neighbor applications	
		if(well_separated(L, R, s)){
			// std::cout << (L->is_leaf() && R->is_leaf()) << std::endl; 
			// std::cout << "well separated true" << std::endl; 
			// std::cout << "prev size of left interactions " << (L->Interactions()).size() << std::endl; 
			// std::cout << "here1" << std::endl; 
			if(L->is_leaf()){
				L -> add_interaction(R);
			}
			// std::cout << "here2" << std::endl; 
			// std::cout << "new size of left interactions " << (L->Interactions()).size() << std::endl; 
			if(R->is_leaf()){
				R -> add_interaction(L);
			}
			// std::cout << "here3" << std::endl; 
		} else{
			size_t n = L->size() + R->size();
			// std::cout << "left size and box " << L->size() << std::endl;
			// print_point(L->Box().first);
			// print_point(L->Box().second);
			// std::cout << "right size and box " << R->size() << std::endl; 
			// print_point(R->Box().first);
			// print_point(R->Box().second);
			double m_L = L->Max_dim();
			// std::cout << "left max dim " << m_L << std::endl; 
			double m_R = R->Max_dim();
			// std::cout << "here4" << std::endl; 
			// std::cout << "right max dim " << m_R << std::endl; 
			if(m_L > m_R){
	  			wsrChildren(R, L->Left(), s);
	  			wsrChildren(R, L->Right(), s);
			} else{
	  			wsrChildren(L, R->Left(), s);
	  			wsrChildren(L, R->Right(), s);
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
    	node* r = build_recursive(parlay::make_slice(P), parlay::make_slice(Tmp));
    	t.next("build");
    	return tree_ptr(r);
  	}


	static void print_point(point p){
		int d = p.dimension();
		std::cout << "Point: ";
		for(int j=0; j<d; j++){
			std::cout << p[j] << ", ";
		}
		std::cout << "\n";
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

  	static node* build_recursive(slice_t points, slice_t Tmp){
  		if(points.size() == 0){
  			std::cout << "ERROR: passed in slice of size 0 when building tree" << std::endl; 
  			abort(); 
  		}
		if(points.size() == 1){
  			return node::new_leaf(points);
  		}else if((points.size()==2) && are_equal(points[0]->pt, points[1]->pt)){
  			node* L = build_recursive(points.cut(0,1), points.cut(0, 1));
  			node* R = build_recursive(points.cut(0,1), points.cut(0, 1));
  			return node::new_node(L, R);
  		}else{
	  		//identify the box around the points and its largest dimension
	  		//TODO see if box calculations can be made more efficient by just splitting
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
	  		// std::cout << d << std::endl; 
	  		double splitpoint = (b.first[d] + b.second[d])/2.0;
	  		// std::cout << splitpoint << std::endl; 
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
	  			// std::cout << nl << std::endl; 
	  		}
	  		// std::cout << "New Round" << std::endl; 
	  		// std::cout << "Left Points" << std::endl; 
	  		// for(size_t i=0; i<nl; i++){
	  		// 	print_point(Tmp[i]->pt);
	  		// }
	  		// std::cout << "Right Points" << std::endl; 
	  		// for(size_t i=nl; i<n; i++){
	  		// 	print_point(Tmp[i]->pt);
	  		// }
	  		// for(size_t i=0; i<n; i++){
	  		// 	print_point(Tmp[i]->pt);
	  		// }
	  		//create left and right children
  			node *L, *R; 
	  		parlay::par_do_if(n > 1000, 
	  			[&] () {L = build_recursive(Tmp.cut(0, nl), points.cut(0, nl));},
	  			[&] () {R = build_recursive(Tmp.cut(nl, n), points.cut(nl, n));}
	  		);
	  		//create parent node
	  		return node::new_node(L, R);
	  	}
  	}

}; //ends CKtree structure

	  		// if(nl==0){
	  		// 	node* R = build_recursive(points.cut(0, nl), Tmp.cut(0, nl));
	  		// }else if(nl==n){
	  		// 	node* L = build_recursive(points.cut(nl, n), Tmp.cut(nl, n)); 
	  		// }else{
	  		// 	node *L, *R; 
		  	// 	parlay::par_do_if(n > 1000, 
		  	// 		[&] () {L = build_recursive(points.cut(0, nl), Tmp.cut(0, nl));},
		  	// 		[&] () {R = build_recursive(points.cut(nl, n), Tmp.cut(nl, n));}
		  	// 	);
	  		// }

