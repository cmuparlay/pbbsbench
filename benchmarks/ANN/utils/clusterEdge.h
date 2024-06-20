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
#include <queue>

/* 
A function which returns a pair of random distinct indices from the active_indices sequence 
*/
std::pair<size_t, size_t> select_two_random(parlay::sequence<size_t>& active_indices,
	parlay::random& rnd) {
	size_t first_index = rnd.ith_rand(0) % active_indices.size(); 
	size_t second_index_unshifted = rnd.ith_rand(1) % (active_indices.size() - 1); // the second index cannot be the last index in active_indices
	size_t second_index = (second_index_unshifted < first_index) ?
	second_index_unshifted : (second_index_unshifted + 1); // if the second index is equal to or greater than the first, it's shifted up by one.
	// is there a reason this isn't a != comparison?

	return {active_indices[first_index], active_indices[second_index]};
}

struct DisjointSet{
	parlay::sequence<int> parent;
	parlay::sequence<int> rank;
	size_t N; 

	DisjointSet(size_t size){
		N = size;
		parent = parlay::sequence<int>(N);
		rank = parlay::sequence<int>(N);
		parlay::parallel_for(0, N, [&] (size_t i) {
			parent[i]=i;
			rank[i] = 0;
		});		
	}

	void _union(int x, int y){
		int xroot = parent[x];
		int yroot = parent[y];
		int xrank = rank[x];
		int yrank = rank[y];
		if(xroot == yroot)
			return;
		else if(xrank < yrank)
			parent[xroot] = yroot;
		else{
			parent[yroot] = xroot;
			if(xrank == yrank)
				rank[xroot] = rank[xroot] + 1;
		}
	}

	int find(int x){
		if(parent[x] != x)
			parent[x] = find(parent[x]);
		return parent[x];
	}

	void flatten(){
		for(int i=0; i<N; i++) find(i);
	}

	bool is_full(){
		flatten();
		parlay::sequence<bool> truthvals(N);
		parlay::parallel_for(0, N, [&] (size_t i){
			truthvals[i] = (parent[i]==parent[0]);
		});
		auto ff = [&] (bool a) {return not a;};
		auto filtered = parlay::filter(truthvals, ff);
		if(filtered.size()==0) return true;
		return false;
	}

};


template<typename T>
struct cluster{
	unsigned d; // dimension of the vectors
	bool mips; // whether comparison is mips as opposed to euclidean
	using tvec_point = Tvec_point<T>;
	using edge = std::pair<int, int>;
	using labelled_edge = std::pair<edge, float>;

	cluster(unsigned dim, bool m): d(dim), mips(m) {}

	float Distance(T* p, T* q, unsigned d){
		if(mips) return mips_distance(p, q, d);
		else return distance(p, q, d);
	}

	/* 
	inserts each edge after checking for duplicates

	@param v: a sequence of tvec_points which are the vertices of the graph
	@param edges: a sequence of edges to be inserted into the graph
	*/
	void process_edges(parlay::sequence<tvec_point*> &v, parlay::sequence<edge> edges){
		int maxDeg = v[1]->out_nbh.begin() - v[0]->out_nbh.begin(); // max degree of a vertex (although I feel like this approach seems a bit fragile when the neighborhoods are not arrays)
		auto grouped = parlay::group_by_key(edges); // group edges by their first index
		for(auto pair : grouped){ // could this be parallelized?
			auto [index, candidates] = pair;
			for(auto c : candidates){
				if(size_of(v[index]->out_nbh) < maxDeg){
					add_nbh(c, v[index]);
				}else{
					remove_edge_duplicates(v[index]);
					add_nbh(c, v[index]);
				}
			}
		}
	}

	void remove_edge_duplicates(tvec_point* p){
		parlay::sequence<int> points;
		for(int i=0; i<size_of(p->out_nbh); i++){
			points.push_back(p->out_nbh[i]);
		}
		auto np = parlay::remove_duplicates(points);
		add_out_nbh(np, p);
	}

	/* 
	bin(N) - bin(N - i)
	*/
	int generate_index(int N, int i){
		return (N*(N-1) - (N-i)*(N-i-1))/2;
	}
	
	/* 
	Builds an approximation of the MSTk of the graph using Kruskal's algorithm

	parameters dim and K are just to interface with the cluster tree code 

	@param v: a sequence of tvec_points which are the vertices of the graph
	@param active_indices: a sequence of indices of the active vertices in v (on which the MSTk is to be built)
	@param dim: the dimension of the vectors
	@param K: the maximum degree of the vertices in the MSTk
	*/
	void MSTk(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices, 
		unsigned dim, int K){
		//preprocessing for Kruskal's
		int N = active_indices.size();
		DisjointSet *disjset = new DisjointSet(N);
		size_t m = 10;
		auto less = [&] (labelled_edge a, labelled_edge b) {return a.second < b.second;};
		parlay::sequence<parlay::sequence<labelled_edge>> pre_labelled(N);
		parlay::parallel_for(0, N, [&] (size_t i){ // for index i of active_indices
			std::priority_queue<labelled_edge, std::vector<labelled_edge>, decltype(less)> Q(less);
			for(int j=0; j<N; j++){ // for every other index j of active_indices
				if(j!=i){
					float dist_ij = Distance(v[active_indices[i]]->coordinates.begin(), v[active_indices[j]]->coordinates.begin(), dim);
					if(Q.size() >= m){
						float topdist = Q.top().second;
						if(dist_ij < topdist){
							labelled_edge e;
							if(i<j) e = std::make_pair(std::make_pair(i,j), dist_ij);
							else e = std::make_pair(std::make_pair(j, i), dist_ij);
							Q.pop();
							Q.push(e);
						}
					}else{
						labelled_edge e;
						if(i<j) e = std::make_pair(std::make_pair(i,j), dist_ij);
						else e = std::make_pair(std::make_pair(j, i), dist_ij);
						Q.push(e);
					}
				}
			}
			parlay::sequence<labelled_edge> edges(m);
			for(int j=0; j<m; j++){edges[j] = Q.top(); Q.pop();}
			pre_labelled[i] = edges;
		});
		auto flat_edges = parlay::flatten(pre_labelled);
		// std::cout << flat_edges.size() << std::endl;
		auto less_dup = [&] (labelled_edge a, labelled_edge b){
			auto dist_a = a.second;
			auto dist_b = b.second;
			if(dist_a == dist_b){
				int i_a = a.first.first;
				int j_a = a.first.second;
				int i_b = b.first.first;
				int j_b = b.first.second;
				if((i_a==i_b) && (j_a==j_b)){
					return true;
				} else{
					if(i_a != i_b) return i_a < i_b;
					else return j_a < j_b;
				}
			}else return (dist_a < dist_b);
		};
		auto labelled_edges = parlay::remove_duplicates_ordered(flat_edges, less_dup);
		// parlay::sort_inplace(labelled_edges, less);
		auto degrees = parlay::tabulate(active_indices.size(), [&] (size_t i) {return 0;});
		parlay::sequence<edge> MST_edges = parlay::sequence<edge>();
		//modified Kruskal's algorithm
		for(int i=0; i<labelled_edges.size(); i++){
			labelled_edge e_l = labelled_edges[i];
			edge e = e_l.first;
			if((disjset->find(e.first) != disjset->find(e.second)) && (degrees[e.first]<K) && (degrees[e.second]<K)){
				MST_edges.push_back(std::make_pair(active_indices[e.first], active_indices[e.second]));
				MST_edges.push_back(std::make_pair(active_indices[e.second], active_indices[e.first]));
				degrees[e.first] += 1;
				degrees[e.second] += 1;
				disjset->_union(e.first, e.second);
			}
			if(i%N==0){
				if(disjset->is_full()){
					break;
				}
			}
		}
		delete disjset;
		process_edges(v, MST_edges);
	}

	bool tvec_equal(tvec_point* a, tvec_point* b, unsigned d){
		for(int i=0; i<d; i++){
			if(a->coordinates[i] != b->coordinates[i]){
				return false;
			}
		}
		return true;
	}

	void recurse(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices,
		parlay::random& rnd, size_t cluster_size, 
		unsigned dim, int K, tvec_point* first, tvec_point* second){
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

		if(closer_first.size() == 1) {
			random_clustering(v, active_indices, right_rnd, cluster_size, dim, K);
		}
		else if(closer_second.size() == 1){
			random_clustering(v, active_indices, left_rnd, cluster_size, dim, K);
		}
		else{
			parlay::par_do(
				[&] () {random_clustering(v, closer_first, left_rnd, cluster_size, dim, K);}, 
				[&] () {random_clustering(v, closer_second, right_rnd, cluster_size, dim, K);}
			);
		}
	}

	void random_clustering(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices,
		parlay::random& rnd, size_t cluster_size, unsigned dim, int K){
		if(active_indices.size() < cluster_size) MSTk(v, active_indices, dim, K);
		else{
			auto [f, s] = select_two_random(active_indices, rnd);
    		tvec_point* first = v[f];
    		tvec_point* second = v[s];

			if(tvec_equal(first, second, dim)){
				// std::cout << "Equal points selected, splitting evenly" << std::endl;
				parlay::sequence<size_t> closer_first;
				parlay::sequence<size_t> closer_second;
				for(int i=0; i<active_indices.size(); i++){
					if(i<active_indices.size()/2) closer_first.push_back(active_indices[i]);
					else closer_second.push_back(active_indices[i]);
				}
				auto left_rnd = rnd.fork(0);
				auto right_rnd = rnd.fork(1);
				parlay::par_do(
					[&] () {random_clustering(v, closer_first, left_rnd, cluster_size, dim, K);}, 
					[&] () {random_clustering(v, closer_second, right_rnd, cluster_size, dim, K);}
				);
			} else{
				recurse(v, active_indices, rnd, cluster_size, dim, K, first, second);
			}
		}
	}

	void random_clustering_wrapper(parlay::sequence<tvec_point*> &v, size_t cluster_size, 
		unsigned dim, int K){
		std::random_device rd;    
  		std::mt19937 rng(rd());   
  		std::uniform_int_distribution<int> uni(0,v.size()); 
    	parlay::random rnd(uni(rng));
    	auto active_indices = parlay::tabulate(v.size(), [&] (size_t i) { return i; });
    	random_clustering(v, active_indices, rnd, cluster_size, dim, K);
	}

	void multiple_clustertrees(parlay::sequence<tvec_point*> &v, size_t cluster_size, int num_clusters,
		unsigned dim, int K, int bound = 0){
		for(int i=0; i<num_clusters; i++){
			random_clustering_wrapper(v, cluster_size, dim, K);
		}
	}
};
