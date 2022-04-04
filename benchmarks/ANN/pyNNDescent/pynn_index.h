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

template<typename T>
struct pyNN_index{
	int K;
	unsigned d;
	using tvec_point = Tvec_point<T>;
	using fine_sequence = parlay::sequence<int>;
	// using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
	using edge = std::pair<int, int>;
	using labelled_edge = std::pair<edge, float>;
	using pid = std::pair<int, float>;

	pyNN_index(int md, unsigned dim) : K(md), d(dim) {}

	static parlay::sequence<edge> naive_neighbors(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t>& active_indices,
		unsigned dim, int maxK){
		size_t n = active_indices.size();
		int k = std::min(maxK, (int) (n-1));
		parlay::sequence<edge> edges(k*n);
		parlay::parallel_for(0, n, [&] (size_t i){
			parlay::sequence<labelled_edge> labelled_edges(n-1);
			parlay::parallel_for(0, n, [&] (size_t j){
				float dist = distance(v[active_indices[i]]->coordinates.begin(), v[active_indices[j]]->coordinates.begin(), dim);
				labelled_edge e = std::make_pair(std::make_pair(active_indices[i], active_indices[j]), dist);
				if(j<i) labelled_edges[j] = e;
				else if(j>i) labelled_edges[j-1] = e;
			});
			auto less = [&] (labelled_edge a, labelled_edge b) {return a.second < b.second;};
			auto sorted_edges = parlay::sort(labelled_edges);
			for(int j=0; j<k; j++) {
				edges[k*i+j] = sorted_edges[j].first;
			}
		});
		return edges;
	}
	
	static parlay::sequence<edge> test(parlay::sequence<tvec_point*> &v, parlay::sequence<size_t>& active_indices,
		unsigned dim, int maxK){
		edge e = std::make_pair(1,2);
		parlay::sequence<edge> test(1);
		test[0] = e;
		return test;
	}

	void build_index(parlay::sequence<tvec_point*> &v, size_t cluster_size){
		clear(v);
		cluster<T> C(d);
		C.multiple_clustertrees(v, cluster_size, 20, &naive_neighbors, d, K);
	}

};

