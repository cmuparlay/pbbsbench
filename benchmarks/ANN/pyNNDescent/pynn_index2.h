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
#include "../utils/clusterPID.h"  
#include "../utils/union.h"  
#include <random>
#include <set>
#include <queue>
#include <math.h>

extern bool report_stats;

template<typename T>
struct pyNN_index{
	int K;
	unsigned d;
	float delta;
	using tvec_point = Tvec_point<T>;
	using fine_sequence = parlay::sequence<int>;
	using edge = std::pair<int, int>;
	using labelled_edge = std::pair<edge, float>;
	using pid = std::pair<int, float>;
	using special_edge = std::pair<int, pid>; //this format makes it easier to use the group_by functions

	parlay::sequence<pid> intermediate_edges;
	parlay::sequence<size_t> intermediate_sizes;

	pyNN_index(int md, unsigned dim, float Delta) : K(md), d(dim), delta(Delta) {}

	parlay::sequence<int> nn_descent(parlay::sequence<tvec_point*> &v, parlay::sequence<int> &changed, 
		int batch_size = 10000){
		size_t n = v.size();
		auto new_changed = parlay::sequence<int>(v.size(), 0);
		//find the edges in the undirected graph
		auto ugraph = truncated_uedges(v, changed);
		std::cout << "generated candidate edges" << std::endl;
		std::pair<int, parlay::sequence<int>> *begin;
		std::pair<int, parlay::sequence<int>> *end = ugraph.begin();
		int counter = 0;
		while(end != ugraph.end()){
			// std::cout << "Round " << counter << " of chunked index" << std::endl;
			counter++;
			begin = end;
			int remaining = ugraph.end() - end;
			end += std::min(remaining, batch_size);
			nn_descent_chunk(v, changed, new_changed, begin, end);
			// std::cout << "Remaining amount " << remaining << std::endl;
		}
		// auto grouped_labelled = parlay::tabulate(ugraph.size(), [&] (size_t i){
		// 	parlay::sequence<special_edge> edges = parlay::sequence<special_edge>();
		// 	for(const int& j: ugraph[i].second){
		// 		for(const int& k: ugraph[i].second){
		// 			if(j!=k){
		// 				float dist = distance(v[j]->coordinates.begin(), v[k]->coordinates.begin(), d);
		// 				edges.push_back(std::make_pair(j, std::make_pair(k, dist)));
		// 				edges.push_back(std::make_pair(k, std::make_pair(j, dist)));
		// 			}
		// 		}
		// 	}
		// 	return edges;
		// });
		// auto grouped_by = parlay::group_by_key(parlay::flatten(grouped_labelled));
		// // update edges of each vertex based on the candidates in grouped_by
		// parlay::parallel_for(0, grouped_by.size(), [&] (size_t i){
		// 	auto [index, candidates] = grouped_by[i];
		// 	for(size_t j=0; j<size_of(v[index]->out_nbh); j++){
		// 		int indexj = v[index]->out_nbh[j];
		// 		float dist = distance(v[index]->coordinates.begin(), v[indexj]->coordinates.begin(), d);
		// 		candidates.push_back(std::make_pair(indexj, dist));
		// 	}
		// 	auto less = [&] (pid a, pid b) {return a.second < b.second;};
		// 	auto sorted_candidates = parlay::remove_duplicates_ordered(candidates, less);
		// 	int k = std::min(K, (int) sorted_candidates.size());
		// 	auto new_nbh = parlay::tabulate(k, [&] (size_t j) {return sorted_candidates[j].first;});
		// 	for(const int& j: new_nbh){
		// 		bool found = false;
		// 		for(size_t k=0; k<size_of(v[index]->out_nbh); k++){
		// 			if(v[index]->out_nbh[k] == j){
		// 				found = true;
		// 				break;
		// 			}
		// 		}
		// 		if(found == false){
		// 			new_changed[index]=1;
		// 			add_out_nbh(new_nbh, v[index]);
		// 			break;
		// 		}
		// 	}
		// });
		return new_changed;	
	}

	void nn_descent_chunk(parlay::sequence<tvec_point*> &v, parlay::sequence<int> &changed, 
		parlay::sequence<int> &new_changed, std::pair<int, parlay::sequence<int>> *begin, 
		std::pair<int, parlay::sequence<int>> *end){
		// tabulate all-pairs distances for the undirected chunk
		size_t stride = end - begin;
		auto less = [&] (pid a, pid b) {return a.second < b.second;};
		// std::cout << "Tabulating distances" << std::endl;
		auto grouped_labelled = parlay::tabulate(stride, [&] (size_t i){
			parlay::sequence<special_edge> edges = parlay::sequence<special_edge>();
			for(const int& j: (begin+i)->second){
				for(const int& k: (begin+i)->second){
					if(j!=k){
						float dist = distance(v[j]->coordinates.begin(), v[k]->coordinates.begin(), d);
						edges.push_back(std::make_pair(j, std::make_pair(k, dist)));
						edges.push_back(std::make_pair(k, std::make_pair(j, dist)));
					}
				}
			}
			return edges;
		});
		// std::cout << "Doing group_by" << std::endl;
		auto grouped_by = parlay::group_by_key(parlay::flatten(grouped_labelled));
		// update edges of each vertex based on the candidates in grouped_by
		// std::cout << "Sorting and updating" << std::endl;
		parlay::parallel_for(0, grouped_by.size(), [&] (size_t i){
			auto [index, candidates] = grouped_by[i];
			size_t indexS = (size_t) index;

			auto sorted_candidates = parlay::sort(candidates, less);

			if(intermediate_sizes[index]>K){
				std::cout << "ERROR: reported size too large" << std::endl;
				abort();
			}
			
			auto merged_candidates = seq_union_pointers(sorted_candidates.begin(), sorted_candidates.end(), 
				intermediate_edges.begin()+indexS*K, intermediate_edges.begin()+indexS*K+intermediate_sizes[index]);


			size_t new_size = std::min(sorted_candidates.size(), (size_t) K);

			for(int j=0; j<new_size; j++){
				bool found = false;
				int cur_index = merged_candidates[j].first;
				for(int k=0; k<intermediate_sizes[indexS]; k++){
					if(intermediate_edges[indexS*K+k].first == cur_index){
						found = true;
						break;
					}
				}
				if(found == false){
					new_changed[index]=1;
					
					intermediate_sizes[indexS] = new_size;
					for(int k=0; k<new_size; k++){
						intermediate_edges[indexS*K+k] = merged_candidates[k];
					}
					break;
				}
			}
		
		});
	}

	int nn_descent_wrapper(parlay::sequence<tvec_point*> &v){
		size_t n = v.size();
		auto changed = parlay::tabulate(n, [&] (size_t i) {return 1;});
		int rounds = 0;
		while(parlay::reduce(changed) >= delta*n){
			std::cout << "Round " << rounds << std::endl; 
			auto new_changed = nn_descent(v, changed);
			changed = new_changed;
			rounds++;
			std::cout << parlay::reduce(changed) << " elements changed" << std::endl; 
		}
		std::cout << "descent converged in " << rounds << " rounds" << std::endl; 
		return rounds;
	}

	parlay::sequence<std::pair<int, parlay::sequence<int>>> truncated_uedges(parlay::sequence<tvec_point*> &v, 
		parlay::sequence<int> &changed){
			
		auto w = parlay::pack_index(changed);
		// std::cout << "computing undirected edges" << std::endl;
		auto grouped_edges = parlay::tabulate(w.size(), [&] (size_t i){
			size_t index = w[i];
			size_t m = intermediate_sizes[index];
			parlay::sequence<special_edge> undirected_edges(2*m);
			for(int j=0; j<m; j++){
				int current_index = intermediate_edges[index*K+j].first;
				float current_dist = intermediate_edges[index*K+j].second;
				undirected_edges[2*j] = std::make_pair(current_index, std::make_pair((int) index, current_dist));
				undirected_edges[2*j+1] = std::make_pair((int) index, std::make_pair(current_index, current_dist));
			}
			return undirected_edges;
		});
		// std::cout << "grouping undirected edges" << std::endl;
		auto ugraph = parlay::group_by_key(parlay::flatten(grouped_edges));
		// std::cout << "trimming edges" << std::endl;
		//trim the undirected edges of each graph
		parlay::sequence<std::pair<int, parlay::sequence<int>>> undgraph(ugraph.size());
		parlay::parallel_for(0, ugraph.size(), [&] (size_t i){
			auto [index, candidates] = ugraph[i];
			auto less = [&] (pid a, pid b) {return a.second < b.second;};
			auto sorted_candidates = parlay::remove_duplicates_ordered(candidates, less);
			int k = std::min((int) sorted_candidates.size(), K*2);
			undgraph[i].first = index;
			undgraph[i].second = parlay::tabulate(k, [&] (size_t j) {return sorted_candidates[j].first;});
		}); 
		return undgraph;
	}


	void undirect_and_prune(parlay::sequence<tvec_point*> &v, double alpha){
		auto grouped_edges = parlay::tabulate(v.size(), [&] (size_t i){
			auto uedges = parlay::tabulate(size_of(v[i]->out_nbh), [&] (size_t j) {
				return std::make_pair(v[i]->out_nbh[j], (int) i);});
			return uedges;
		});
		auto reverse_graph = parlay::group_by_key(parlay::flatten(grouped_edges));
		
		auto less = [&] (pid a, pid b) {return a.second < b.second;};
		// calculate distances, prune triangles, truncate
		parlay::parallel_for(0, reverse_graph.size(), [&] (size_t i){
			auto [index, reverse_candidates] = reverse_graph[i];
			auto reverse_distances = parlay::tabulate(reverse_candidates.size(), [&] (size_t j){
				float dist = distance(v[index]->coordinates.begin(), v[reverse_candidates[j]]->coordinates.begin(), d);
				return std::make_pair(reverse_candidates[j], dist);
			});
			auto distances = parlay::tabulate(size_of(v[index]->out_nbh), [&] (size_t j){
				float dist = distance(v[index]->coordinates.begin(), v[v[index]->out_nbh[j]]->coordinates.begin(), d);
				return std::make_pair(v[index]->out_nbh[j], dist);
			});

			auto sorted_reverse = parlay::sort(reverse_distances, less);
			auto sorted = parlay::sort(distances, less);
			auto candidates = seq_union(sorted_reverse, sorted);
			parlay::sequence<int> new_out = parlay::sequence<int>();
			for(const pid& j : candidates){
				if(new_out.size() == 2*K) break;
				else if(new_out.size() == 0) new_out.push_back(j.first);
				else{
					float dist_p = j.second;
					bool add = true;
					for(const int& k : new_out){
						float dist = distance(v[j.first]->coordinates.begin(), v[k]->coordinates.begin(), d);
						if(dist_p > alpha*dist) {add = false; break;}
					}
					if(add) new_out.push_back(j.first);
				}
			}
			add_out_nbh(new_out, v[index]);
		});

	}

	void assign_edges(parlay::sequence<tvec_point*> &v){
		parlay::parallel_for(0, v.size(), [&] (size_t i){
			for(size_t j=0; j<intermediate_sizes[i]; j++){
				int out_nbh = intermediate_edges[i*K+j].first;
				add_nbh(out_nbh, v[i]);
			}
		});
	}


	void build_index(parlay::sequence<tvec_point*> &v, size_t cluster_size, int num_clusters, double alpha){
		clear(v);
		intermediate_edges = parlay::sequence<pid>(v.size()*K);
        intermediate_sizes = parlay::sequence<size_t>(v.size(), 0);
		clusterPID<T> C(d);
		C.multiple_clustertrees(v, cluster_size, num_clusters, d, K, intermediate_edges, intermediate_sizes);
		std::cout << "finished clustering" << std::endl; 
		nn_descent_wrapper(v);
		assign_edges(v);
		undirect_and_prune(v, alpha);
		// std::cout << "prepared final graph" << std::endl;
	}

};

