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
#include <set>
#include <unordered_set>
#include <random>
#include "common/geometry.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "../utils/types.h"

template<typename T>
struct cluster_tree{
	using tvec_point = Tvec_point<T>;
	struct node{
		node* R;
		node* L;
		bool is_leaf; 
		int center;
		parlay::sequence<size_t> elts;
		size_t n; 

		//leaf node
		node(parlay::sequence<size_t> Elts, int ctr) : elts(Elts), center(ctr) {is_leaf = true; n = Elts.size();}

		//internal node
		node(node* right, node* left, int ctr): R(right), L(left), center(ctr) {is_leaf = false; n = R->n + L->n;}

		static node* new_leaf(parlay::sequence<size_t> elts, int ctr) {
	      node* r = alloc_node();
	      new (r) node(elts, ctr);
	      return r;
	    }

	    static node* new_node(node* L, node* R, int ctr) {
	      node* nd = alloc_node();
	      new (nd) node(L, R, ctr);
	      return nd;
	    }
	    
	    ~node() {
	      // need to collect in parallel
	      parlay::par_do_if(n > 1000,
				[&] () { delete_tree(L);},
				[&] () { delete_tree(R);});
	    }

	    static void delete_tree(node* N) {
	    	if (N != nullptr) {
				N->~node();
				node::free_node(N);
	      	}	
	    }

		private:
			static parlay::type_allocator<node> node_allocator;
		    static node* alloc_node() {
		      return node_allocator.alloc();}
		    static void free_node(node* N) {
		      node_allocator.free(N);}
	};

	static std::pair<size_t, size_t> select_two_random(parlay::sequence<size_t>& active_indices,
      	parlay::random& rnd) {
    	size_t first_index = rnd.ith_rand(0) % active_indices.size(); 
    	size_t second_index_unshifted = rnd.ith_rand(1) % (active_indices.size()-1);
    	size_t second_index = (second_index_unshifted < first_index) ?
      	second_index_unshifted : (second_index_unshifted + 1);
    	return {active_indices[first_index], active_indices[second_index]};
  	}

	static node* random_clustering(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices,
		parlay::random& rnd, size_t cluster_size, int center, unsigned d){
		if(active_indices.size() <= cluster_size) return node::new_leaf(active_indices, center);
		else{
			auto [f, s] = select_two_random(active_indices, rnd);
    		tvec_point* first = v[f];
    		tvec_point* second = v[s];

    		// Split points based on which of the two points are closer.
		    auto closer_first = parlay::filter(parlay::make_slice(active_indices), [&] (size_t ind) {
		      tvec_point* p = v[ind];
		      float dist_first = distance(p->coordinates.begin(), first->coordinates.begin(), d);
		      float dist_second = distance(p->coordinates.begin(), second->coordinates.begin(), d);
		      return dist_first <= dist_second;

		    });

		    auto closer_second = parlay::filter(parlay::make_slice(active_indices), [&] (size_t ind) {
		      tvec_point* p = v[ind];
		      float dist_first = distance(p->coordinates.begin(), first->coordinates.begin(), d);
		      float dist_second = distance(p->coordinates.begin(), second->coordinates.begin(), d);
		      return dist_second < dist_first;
		    });

		    auto left_rnd = rnd.fork(0);
		    auto right_rnd = rnd.fork(1);

			node* L;
			node* R;
			parlay::par_do(
				[&] () {L = random_clustering(v, closer_first, left_rnd, cluster_size, f, d);}, 
				[&] () {R = random_clustering(v, closer_second, right_rnd, cluster_size, s, d);}
			);
			return node::new_node(L, R, center);
		}
	}

	static node* random_clustering_wrapper(parlay::sequence<tvec_point*> &v, size_t cluster_size, unsigned d){
		std::random_device rd;    
  		std::mt19937 rng(rd());   
  		std::uniform_int_distribution<int> uni(0,v.size()); 
    	parlay::random rnd(uni(rng));
    	// std::cout << "tabulating" << std::endl; 
    	auto active_indices = parlay::tabulate(v.size(), [&] (size_t i) { return i; });
    	//give a center of 0 to the root because it doesn't matter
    	// std::cout << "clustering" << std::endl;
    	return random_clustering(v, active_indices, rnd, cluster_size, 0, d); 
	}

	static std::pair<size_t, float> search(node* N, tvec_point* p, parlay::sequence<tvec_point*> &v, unsigned d){
		if(N->is_leaf){
			auto candidates = parlay::tabulate(N->elts.size(), [&] (size_t i) {
				float dist = distance(p->coordinates.begin(), v[N->elts[i]]->coordinates.begin(), d);
				return std::make_pair(N->elts[i], dist);
			});
			auto comp = [&] (std::pair<size_t, float> a, std::pair<size_t, float> b) {return a.second < b.second;};
			return *(parlay::min_element(candidates, comp));
		} else{
			size_t ctr_l = (N->L)->center;
			size_t ctr_r = (N->R)->center;
			float d_l = distance(p->coordinates.begin(), v[ctr_l]->coordinates.begin(), d);
			float d_r = distance(p->coordinates.begin(), v[ctr_r]->coordinates.begin(), d);
			if(d_l<d_r) return search(N->L, p, v, d);
			else return search(N->R, p, v, d);
		}
	}

};

