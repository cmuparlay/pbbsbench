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
#include "common/geometry.h"
#include "indexTools.h"
#include <set>

template<typename T>
void graph_stats(parlay::sequence<Tvec_point<T>*> &v){
	auto od = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return size_of(v[i]->out_nbh);});
  	size_t j = parlay::max_element(od) - od.begin();
  	int maxDegree = od[j];
  	size_t k = parlay::min_element(od) - od.begin();
  	size_t sum1 = parlay::reduce(od);
  	std::cout << "Average graph out-degree = " << sum1/((double) v.size()) << std::endl;
  	std::cout << "Max out-degree: " << maxDegree << ", Min out-degree: " << od[k] << std::endl;  
}

template<typename T>
void query_stats(parlay::sequence<Tvec_point<T>*> &q){
	//stats on visited lists
	auto visited_stats = parlay::tabulate(q.size(), [&] (size_t i) {return q[i]->visited;});
	parlay::sort_inplace(visited_stats);
	int avg_visited = (int) parlay::reduce(visited_stats)/((double) q.size());
	size_t tail_index = .99*((float) q.size());
	int tail_visited = visited_stats[tail_index];
	std::cout << "Average num visited: " << avg_visited << ", 99th percentile num visited: " << tail_visited << std::endl;
	//stats on distance comparisons
	auto dist_stats = parlay::tabulate(q.size(), [&] (size_t i) {return q[i]->dist_calls;});
	parlay::sort_inplace(dist_stats);
	int avg_dist = (int) parlay::reduce(dist_stats)/((double) q.size());
	int tail_dist = dist_stats[tail_index];
	std::cout << "Average dist cmps: " << avg_dist << ", 99th percentile dist cmps: " << tail_dist << std::endl;
	std::cout << std::endl;
}