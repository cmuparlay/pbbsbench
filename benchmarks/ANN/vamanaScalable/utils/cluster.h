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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "common/geometry.h"
#include "indexTools.h"
#include <random>
#include <set>
#include <math.h>
#include <functional>

template<typename T>
struct cluster{
	unsigned d; 
	using tvec_point = Tvec_point<T>;
	using edge = std::pair<int, int>;
	using labelled_edge = std::pair<edge, float>;

	cluster(unsigned dim): d(dim){}

	std::pair<size_t, size_t> select_two_random(parlay::sequence<size_t>& active_indices,
      	parlay::random& rnd) {
    	size_t first_index = rnd.ith_rand(0) % active_indices.size(); 
    	size_t second_index_unshifted = rnd.ith_rand(1) % (active_indices.size()-1);
    	size_t second_index = (second_index_unshifted < first_index) ?
      	second_index_unshifted : (second_index_unshifted + 1);

    	return {active_indices[first_index], active_indices[second_index]};
  	}

	parlay::sequence<edge> random_clustering(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices,
		parlay::random& rnd, size_t cluster_size, 
		function<parlay::sequence<edge>(parlay::sequence<tvec_point*>&, parlay::sequence<size_t>&, unsigned, int)> edge_generator, 
		unsigned dim, int K){
		if(active_indices.size() < cluster_size) return edge_generator(v, active_indices, dim, K);
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

			parlay::sequence<edge> left_edges;
			parlay::sequence<edge> right_edges;
			if(closer_first.size() == 1) {
				return random_clustering(v, closer_second, right_rnd, cluster_size, edge_generator, dim, K);
			}
			else if(closer_second.size() == 1){
				return left_edges = random_clustering(v, closer_first, left_rnd, cluster_size, edge_generator, dim, K);
			}
			else{
				parlay::par_do(
					[&] () {left_edges = random_clustering(v, closer_first, left_rnd, cluster_size, edge_generator, dim, K);}, 
					[&] () {right_edges = random_clustering(v, closer_second, right_rnd, cluster_size, edge_generator, dim, K);}
				);
				parlay::sequence<parlay::sequence<edge>> to_flatten(2);
				to_flatten[0] = left_edges;
				to_flatten[1] = right_edges;
				return parlay::flatten(to_flatten);
			}
		}
	}

	parlay::sequence<edge> random_clustering_wrapper(parlay::sequence<tvec_point*> &v, size_t cluster_size, 
		function<parlay::sequence<edge>(parlay::sequence<tvec_point*>&, parlay::sequence<size_t>&, unsigned, int)> edge_generator, 
		unsigned dim, int K){
		std::random_device rd;    
  		std::mt19937 rng(rd());   
  		std::uniform_int_distribution<int> uni(0,v.size()); 
    	parlay::random rnd(uni(rng));
    	auto active_indices = parlay::tabulate(v.size(), [&] (size_t i) { return i; });
    	return random_clustering(v, active_indices, rnd, cluster_size, edge_generator, dim, K);
	}

	//takes in multiple sequences of edges and assigns them to points in v
	void edge_union(parlay::sequence<tvec_point*> &v, parlay::sequence<parlay::sequence<edge>> &edges){
		size_t n = v.size();
		auto grouped_by = parlay::group_by_key(parlay::flatten(edges));
		parlay::parallel_for(0, n, [&] (size_t i){
			auto [index, candidates] = grouped_by[i];
			auto nbh = parlay::remove_duplicates(candidates);
			add_out_nbh(nbh, v[index]);
		});
	}

	void multiple_clustertrees(parlay::sequence<tvec_point*> &v, size_t cluster_size, int num_clusters,
		function<parlay::sequence<edge>(parlay::sequence<tvec_point*>&, parlay::sequence<size_t>&, unsigned, int)> edge_generator, 
		unsigned dim, int K){
		parlay::sequence<parlay::sequence<edge>> edge_wrapped(num_clusters);
			parlay::parallel_for(0, num_clusters, [&] (size_t i){
				edge_wrapped[i] = random_clustering_wrapper(v, cluster_size, edge_generator, dim, K);
			});
		edge_union(v, edge_wrapped);
	}
};