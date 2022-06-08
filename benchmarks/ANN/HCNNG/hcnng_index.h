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
#include "../utils/cluster.h"
#include <random>
#include <set>
#include <math.h>

extern bool report_stats;

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
struct hcnng_index{
	int maxDeg;
	unsigned d;
	using tvec_point = Tvec_point<T>;
	using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
	using edge = std::pair<int, int>;
	using labelled_edge = std::pair<edge, float>;
	using pid = std::pair<int, float>;

	hcnng_index(int md, unsigned dim) : maxDeg(md), d(dim) {}

	static int generate_index(int N, int i){
		return (N*(N-1) - (N-i)*(N-i-1))/2;
	}

	//parameters dim and K are just to interface with the cluster tree code
	static parlay::sequence<edge> MSTk(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t> &active_indices, 
		unsigned dim, int K){
		//preprocessing for Kruskal's
		int N = active_indices.size();
		DisjointSet *disjset = new DisjointSet(N);
		parlay::sequence<labelled_edge> labelled_edges(N*(N-1)/2);
		parlay::parallel_for(0, N, [&] (size_t i) {
			size_t start = generate_index(N, i);
			parlay::parallel_for(i+1, N, [&] (size_t j) {
				float dist_ij = distance(v[active_indices[i]]->coordinates.begin(), v[active_indices[j]]->coordinates.begin(), dim);
				labelled_edges[start+j-(i+1)] = std::make_pair(std::make_pair(i,j), dist_ij);
			});
			
		});
		auto less = [&] (labelled_edge a, labelled_edge b) {return a.second < b.second;};
		auto sorted_edges = parlay::sort(labelled_edges, less);
		auto degrees = parlay::tabulate(active_indices.size(), [&] (size_t i) {return 0;});
		parlay::sequence<edge> MST_edges = parlay::sequence<edge>();
		//modified Kruskal's algorithm
		for(int i=0; i<sorted_edges.size(); i++){
			labelled_edge e_l = sorted_edges[i];
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
		return MST_edges;
	}


	void build_index(parlay::sequence<tvec_point*> &v, int cluster_rounds, size_t cluster_size){ 
		clear(v); 
		cluster<T> C(d);
		C.multiple_clustertrees(v, cluster_size, cluster_rounds, &MSTk, d, maxDeg);
	}
	
};