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
#include "clusterEdge.h"
#include "union.h"
#include <random>
#include <set>
#include <math.h>
#include <functional>


template<typename T>
struct clusterPID{
	unsigned d; 
	using tvec_point = Tvec_point<T>;
	using edge = std::pair<int, int>;
    using pid = std::pair<int, float>;
    // parlay::sequence<pid> intermediate_edges;
    // parlay::sequence<size_t> intermediate_sizes;

	clusterPID(unsigned dim): d(dim){}

    void naive_neighbors(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t>& active_indices,
		unsigned dim, int maxK, parlay::sequence<pid> &intermediate_edges, 
		parlay::sequence<size_t> &intermediate_sizes){
		size_t n = active_indices.size();
		parlay::parallel_for(0, n, [&] (size_t i){
			parlay::sequence<pid> candidates;
            size_t index = active_indices[i];
			//tabulate all-pairs distances between the elements in the leaf
			for(int j=0; j<i; j++){
                if(j != i){
                    float dist = distance(v[index]->coordinates.begin(), v[active_indices[j]]->coordinates.begin(), dim);
                    pid e = std::make_pair(active_indices[j], dist);
                    candidates.push_back(e);
                }		
			}
			auto less = [&] (pid a, pid b) {return a.second < b.second;};
			auto sorted_edges = parlay::sort(candidates);
            size_t intermediate_size = intermediate_sizes[index];
            auto new_best = seq_union_pointers(candidates.begin(), candidates.end(), 
                intermediate_edges.begin() + (int) maxK*index, 
                intermediate_edges.begin() + (int) maxK*index+intermediate_size);
			size_t max_size = std::min((size_t) maxK, new_best.size());
			for(int j=0; j<max_size; j++){
                intermediate_edges[maxK*index+j] = new_best[j];
            }
            intermediate_sizes[index] = max_size;
		});
	}
	
	bool tvec_equal(tvec_point* a, tvec_point* b, unsigned d){
		for(int i=0; i<d; i++){
			if(a->coordinates[i] != b->coordinates[i]){
				return false;
			}
		}
		return true;
	}

	void random_clustering(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices,
		parlay::random& rnd, size_t cluster_size, unsigned dim, int K, parlay::sequence<pid> &intermediate_edges, 
		parlay::sequence<size_t> &intermediate_sizes){
		if(active_indices.size() < cluster_size)  naive_neighbors(v, active_indices, dim, K, 
			intermediate_edges, intermediate_sizes);
		else{
			auto [f, s] = select_two_random(active_indices, rnd);
    		tvec_point* first = v[f];
    		tvec_point* second = v[s];

			auto left_rnd = rnd.fork(0);
			auto right_rnd = rnd.fork(1);


			if(tvec_equal(first, second, dim)){
				std::cout << "Equal points selected, splitting evenly" << std::endl;
				parlay::sequence<size_t> closer_first;
				parlay::sequence<size_t> closer_second;
				for(int i=0; i<active_indices.size(); i++){
					if(i<active_indices.size()/2) closer_first.push_back(active_indices[i]);
					else closer_second.push_back(active_indices[i]);
				}
				auto left_rnd = rnd.fork(0);
				auto right_rnd = rnd.fork(1);
				parlay::par_do(
					[&] () {random_clustering(v, closer_first, left_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);}, 
					[&] () {random_clustering(v, closer_second, right_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);}
				);
			} else{
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

				if(closer_first.size() == 1) {
					random_clustering(v, closer_second, right_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);
				}
				else if(closer_second.size() == 1){
					random_clustering(v, closer_first, left_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);
				}
				else{
					parlay::par_do(
						[&] () {random_clustering(v, closer_first, left_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);}, 
						[&] () {random_clustering(v, closer_second, right_rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);}
					);
				}
			}
		}
	}

	void random_clustering_wrapper(parlay::sequence<tvec_point*> &v, size_t cluster_size, 
		unsigned dim, int K, parlay::sequence<pid> &intermediate_edges, 
		parlay::sequence<size_t> &intermediate_sizes){
		std::random_device rd;    
  		std::mt19937 rng(rd());   
  		std::uniform_int_distribution<int> uni(0,v.size()); 
    	parlay::random rnd(uni(rng));
    	auto active_indices = parlay::tabulate(v.size(), [&] (size_t i) { return i; });
    	random_clustering(v, active_indices, rnd, cluster_size, dim, K, intermediate_edges, intermediate_sizes);
	}

	void multiple_clustertrees(parlay::sequence<tvec_point*> &v, size_t cluster_size, int num_clusters,
		unsigned dim, int K, parlay::sequence<pid> &intermediate_edges, 
		parlay::sequence<size_t> &intermediate_sizes){
		for(int i=0; i<num_clusters; i++){
			std::cout << "Cluster " << i << std::endl;
			random_clustering_wrapper(v, cluster_size, dim, K, intermediate_edges, intermediate_sizes);
		}
	}
};