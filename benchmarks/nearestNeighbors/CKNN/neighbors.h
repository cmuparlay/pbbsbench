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
#include <math.h> 
#include <queue>
#include <limits>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/hash_table.h"
#include "common/geometry.h"
#include "ck_tree.h"

bool report_stats = true; 
bool check_correctness = false; 
int algorithm_version = 0; //just keep this at 0

template<class vtx, int max_k>
struct k_nearest_neighbors{

	using point = typename vtx::pointT;
	using ck_tree_vtx = ck_tree<vtx>; 
	using node = typename ck_tree_vtx::node; 
	using tree_ptr = typename ck_tree_vtx::tree_ptr;
  	using box = typename ck_tree_vtx::box;
  	using slice_n = decltype(make_slice(parlay::sequence<node*>()));
  	using vtx_dist = std::pair<vtx*, double>;

  	tree_ptr tree; 

  	k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    	tree = ck_tree_vtx::build(V); 
  	}

  	void get_wsr(node* T, double s, int k){
  		ck_tree_vtx::wsr(T, s, k);
  	}

	struct kNN{

		vtx *vertex;  // the vertex for which we are trying to find a NN
		node *leaf; //the leaf containing the vertex 
		int k;
    	vtx *neighbors[max_k];  // the current k nearest neighbors (nearest last)
    	double distances[max_k]; // distance to current k nearest neighbors
    	int interactions_cnt; 

    	kNN() {}


	    // returns the ith smallest element (0 is smallest) up to k-1
	    vtx* operator[] (const int i) { return neighbors[k-i-1]; } 
	    
	    kNN(vtx *p, node* T, int kk) {
		    if (kk > max_k) {
				std::cout << "k too large in kNN" << std::endl;
				abort();
			}
		    k = kk;
		    vertex = p;
		    leaf = T; 
		    // initialize nearest neighbors to point to Null with
		    // distance "infinity".
	        for (int i=0; i<k; i++) {
	        	neighbors[i] = (vtx*) NULL; 
	        	distances[i] = numeric_limits<double>::max();
	        }
	    }
	    

    	slice_n flatten(parlay::sequence<node*> interactions){
    		size_t size = interactions.size();
    		size_t num_leaves = 0;
    		for(size_t i=0; i<size; i++){
    			num_leaves += interactions[i]->size(); 
    		}
    		parlay::sequence<node*> leaves;
    		leaves = parlay::sequence<node*>(num_leaves);
    		slice_n leaves_slice = parlay::make_slice(leaves);
    		size_t last_idx = 0; 
    		for(size_t i=0; i<size; i++){
    			int num_children = interactions[i]->size();
    			flatten_rec(interactions[i], leaves_slice.cut(last_idx, last_idx+num_children));
    			last_idx += num_children; 
    		}
    		return leaves_slice; 
    	}

    	void flatten_rec(node* T, slice_n ints){
    		if(T->is_leaf()){
    			ints[0] = T;
    		} else{
    			node* L = T->Left();
    			node* R = T->Right();
    			int num_leaves_left = L->size();
    			int num_leaves_right = R->size();
    			int num_leaves = num_leaves_left + num_leaves_right;
    			// parlay::par_do_if(num_leaves > 2,
    			// 	[&] () {flatten_rec(L, ints.cut(0, num_leaves_left));},
    			// 	[&] () {flatten_rec(R, ints.cut(num_leaves_left, num_leaves_left+num_leaves_right));}
    			// );
				flatten_rec(L, ints.cut(0, num_leaves_left));
				flatten_rec(R, ints.cut(num_leaves_left, num_leaves_left+num_leaves_right));
    		}
    	}

    	void update_nonrec(slice_n leaves){
    		for(size_t i=0; i<leaves.size(); i++){
    			update_nearest(leaves[i]->Vertex());
    		}
    	}

    	void update_nearest(vtx *other) {  
      		auto dist = (vertex->pt - other->pt).sqLength();
      		if (dist < distances[0]) { 
      			neighbors[0] = other;
      			distances[0] = dist;
	      		for (int i = 1;
	      	     	i < k && distances[i-1] < distances[i];
	      	     	i++) {
	      	  		swap(distances[i-1], distances[i]);
	      	  		swap(neighbors[i-1], neighbors[i]); 
	        	}
      		}
    	}

    	vtx_dist update_rec(slice_n leaves){
    		size_t num_leaves = leaves.size();
    		if (num_leaves==1){
    			auto dist = (vertex->pt - (leaves[0]->Vertex())->pt).sqLength();
    			return std::make_pair(leaves[0]->Vertex(), dist);
    		} else{
    			vtx_dist L, R; 
    			L = update_rec(leaves.cut(0, num_leaves/2));
    			R = update_rec(leaves.cut(num_leaves/2, num_leaves));
    			if(L.second < R.second) return L;
    			else return R; 
    		}
    	}

    	parlay::sequence<node*> get_interactions(){
    		node* current = leaf; 
    		parlay::sequence<node*> ints_per_node = parlay::sequence<node*>();
    		while(current->size() <= k && current->Parent() != nullptr){
    			ints_per_node.push_back(current);
    			current = current->Parent();
    		}
    		parlay::sequence<node*> all_ints = parlay::sequence<node*>();
    		for(size_t i=0; i<ints_per_node.size(); i++){
    			parlay::sequence<node*> interactions = ints_per_node[i]->Interactions();
    			for(size_t j=0; j<interactions.size(); j++){
    				all_ints.push_back(interactions[j]);
    			}
    		}
    		return all_ints;
    	}

    	void k_nearest_neighbors(){
    		parlay::sequence<node*> interactions = get_interactions();
    		slice_n leaves = flatten(interactions);
    		if (report_stats) vertex->counter = leaves.size();
    		update_nonrec(leaves);
    	}

	}; //end KNN struct

  	//brute-force check correctness for approximately 1% of the points
	bool do_check_correct(size_t n){
		float check = (float) rand()/RAND_MAX;
		float threshold = (1/(n*.01));
		if (check < threshold) return true;
		else return false;
	}

	bool are_equal(point p, point q){
		int d = p.dimension();
		for (int j = 0; j<d; j++){
			if(p[j] != q[j]){
				return false;
			}
		}
		return true;
	}

	static void print_point(point p){
		int d = p.dimension();
		std::cout << "Point: ";
		for(int j=0; j<d; j++){
			std::cout << p[j] << ", ";
		}
		std::cout << "\n";
	}

	void update_nearest(vtx* other, vtx* vertex, vtx* neighbors[], double distances[], int k) {  
  		auto dist = (vertex->pt - other->pt).sqLength();
  		if (dist < distances[0]) { 
  			neighbors[k] = other;
  			distances[k] = dist;
      		for (int i = 1;
      	     	i < k && distances[k-i] < distances[k-i+1];
      	     	i++) {
      	  		swap(distances[k-i], distances[k-i+1]);
      	  		swap(neighbors[k-i], neighbors[k-i+1]); 
        	}
  		}
	}

	void check_correct_vtx(vtx* v, parlay::sequence<vtx*> V, int k){
		size_t n = V.size(); 
		if (do_check_correct(n)){
			//collect reported neighbors and distances
			vtx* original_neighbors[k]; 
			double original_distances[k]; 
			for(int i=0; i<k; i++){
				original_neighbors[i] = v->ngh[i];
				original_distances[k] = ((v->pt) - (original_neighbors[i]->pt)).sqLength();
			}
			//initialize vectors for correctness check
			vtx* neighbors[k];  
    		double distances[k];
    		for (int i=0; i<k; i++) {
	        	neighbors[i] = (vtx*) NULL; 
	        	distances[i] = numeric_limits<double>::max();
	        }
			for(size_t i = 0; i < n; i++){
				if(not are_equal(V[i]->pt, v->pt)){ //make sure we don't report the query point as its own nearest neighbor
					update_nearest(V[i], v, neighbors, distances, k);
				}
			} 
			for(int i=0; i<k; i++){
				if(not (original_distances[i] <= distances[i])){
					std::cout<<"ERROR: nearest neighbor not correct for neighbor " << i << "\n";
					std::cout << "Query point: ";
					print_point(v->pt);
					std::cout << "Reported neighbor: ";
					print_point(original_neighbors[i]->pt);
					std::cout << "Reported distance: " << (original_neighbors[i]->pt - v->pt).sqLength() << "\n";
					std::cout << "Actual neighbor: ";
					print_point(neighbors[i]->pt);
					std::cout << "Actual distance: " << distances[i] << "\n";
					abort();
				}
			}
		}		
	}

	void check_correct(parlay::sequence<vtx*> V, int k){
		size_t n = V.size();
		for(size_t i=0; i<n; i++){
			check_correct_vtx(V[i], V, k);
		}
	}

	void k_nearest(vtx* v, node* T, int k){ //, parlay::sequence<vtx*> V){
		kNN nn(v, T, k);
		nn.k_nearest_neighbors();
		for(int i=0; i<k; i++){
			v->ngh[i] = nn.neighbors[i];
		}
	}

}; //end k_nearest_neighbors struct



double s=2.1; 

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

	T.get_wsr(T.tree.get(), s, k);
	t.next("find well-separated realization");  


	auto f = [&] (vtx* p, node* n){
		return T.k_nearest(p, n, k);
	};

	T.tree->map(f);
	t.next("compute nearest neighbors");
	if(check_correctness){
		T.check_correct(v, k);
		t.next("check correctness");
	}

	if (report_stats) {
		auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
      	size_t i = parlay::max_element(s) - s.begin();
      	size_t sum = parlay::reduce(s);
      	std::cout << "max interactions = " << s[i] 
			<< ", average interactions = " << sum/((double) v.size()) << std::endl;
		t.next("stats");
	}
	t.next("delete tree");
  };


}

