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
#include "common/geometry.h"
#include "ck_tree.h"

bool report_stats = true; 
int algorithm_version = 0; //just keep this at 0

template<class vtx, int max_k>
struct k_nearest_neighbors{

	using point = typename vtx::pointT;
	using ck_tree_vtx = ck_tree<vtx>; 
	using node = typename ck_tree_vtx::node; 
	using tree_ptr = typename ck_tree_vtx::tree_ptr;
  	using box = typename ck_tree_vtx::box;
  	using slice_n = decltype(make_slice(parlay::sequence<node*>()));

  	tree_ptr tree; 

  	k_nearest_neighbors(parlay::sequence<vtx*> &V) {
    	tree = ck_tree_vtx::build(V); 
  	}

  	void get_wsr(node* T, double s){
  		ck_tree_vtx::wsr(T, s);
  	}

	struct kNN{

		vtx *vertex;  // the vertex for which we are trying to find a NN
		node *leaf; //the leaf containing the vertex 
		int k;
    	vtx *neighbors[max_k];  // the current k nearest neighbors (nearest last)
    	double distances[max_k]; // distance to current k nearest neighbors

    	kNN() {}

    	vtx* get_nbh() {return neighbors[0];}
    	double get_dist() {return distances[0];}

	    // returns the ith smallest element (0 is smallest) up to k-1
	    vtx* operator[] (const int i) { return neighbors[k-i-1]; } 
	    
	    kNN(vtx *p, node* T, int kk) {
		    if (kk > max_k) {
				std::cout << "k too large in kNN" << std::endl;
				abort();
			}
		    k = kk;
		    // print_point(p->pt);
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
    			flatten_rec(L, ints.cut(0, num_leaves_left));
    			flatten_rec(R, ints.cut(num_leaves_left, num_leaves_left+num_leaves_right));
    		}
    	}

    	void update_nearest(vtx *other) {  
    		// std::cout << "updating nearest" << std::endl; 
    		// print_point(vertex->pt);
    		// print_point(other->pt);
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

    	void k_nearest_neighbors(){
    		parlay::sequence<node*> interactions = leaf->Interactions();
    		// std::cout << "size of interaction list " << interactions.size() << std::endl; 
    		slice_n leaves = flatten(interactions);
    		size_t num_leaves = leaves.size();
    		// std::cout << "number of leaves being searched " << num_leaves << std::endl; 
    		// std::cout << "query point:" << std::endl; 
    		// print_point(vertex->pt);
    		// std::cout << "interactions:" << std::endl; 
    		for(size_t i=0; i<num_leaves; i++){
    			// print_point((leaves[i]->Vertex())->pt);
    			update_nearest(leaves[i]->Vertex());
    		}
    	}

	};

  	//brute-force check correctness for approximately 100 points out of every 10 million
	bool do_check_correct(){
		// float check = (float) rand()/RAND_MAX;
		// if (check < .00001) return true;
		// return false;
		return true; 
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

	void check_correct(vtx* v, vtx* ans, double original_dist, parlay::sequence<vtx*> V){
		// print_point(ans->pt);
		// std::cout << original_dist << std::endl;
		if (do_check_correct()){
			// std::cout << "here1" << std::endl; 
			size_t n = V.size(); 
			vtx* nearest;
			double nearest_dist = numeric_limits<double>::max();
			// std::cout << "here2" << std::endl; 
			for(size_t i = 0; i < n; i++){
				// std::cout << "here3" << std::endl; 
				if(not are_equal(V[i]->pt, v->pt)){ //make sure we don't report the query point as its own nearest neighbor
					// std::cout << "here4" << std::endl; 
					double dist = (v->pt - V[i]->pt).sqLength();
					if(dist < nearest_dist){
						nearest = V[i];
						nearest_dist = dist; 
					}
				}
			} 
			if(not (original_dist <= nearest_dist)){
				std::cout << "Query point: ";
				print_point(v->pt);
				std::cout << "Reported neighbor: ";
				print_point(ans->pt);
				std::cout << "Reported distance: " << (ans->pt - v->pt).sqLength() << "\n";
				std::cout << "Actual neighbor: ";
				print_point(nearest->pt);
				std::cout << "Actual distance: " << nearest_dist << "\n";
				std::cout<<"ERROR: nearest neighbor not correct"<< "\n";
				abort();
			}
		}		
	}

	void k_nearest(vtx* v, node* T, int k, parlay::sequence<vtx*> V){
		kNN nn(v, T, k);
		nn.k_nearest_neighbors();
		vtx* ans = nn.get_nbh();
		double dist = nn.get_dist();
		check_correct(v, ans, dist, V);
	}

	void print_tree(){
		node* root = tree.get();
		node* L = root->Left();
		node* R = root->Right();
		std::cout << "Left children: " << std::endl; 
		print_point(((L->Left())->Vertex())->pt);
		print_point(((L->Right())->Vertex())->pt);
		std::cout << "Right children: " << std::endl; 
		print_point(((R->Left())->Vertex())->pt);
		print_point(((R->Right())->Vertex())->pt);
	}

};



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

	T.get_wsr(T.tree.get(), s);
	t.next("find well-separated realization");

	// auto f = [&] (vtx* p, node* n){
	// 	return T.k_nearest(p, n, k, v);
	// };

	// T.tree->map(f);
	// T.print_tree(); 
	t.next("compute nearest neighbors");

	// if (report_stats) {
	// t.next("stats");
	// }
	// t.next("delete tree");
  };


}

