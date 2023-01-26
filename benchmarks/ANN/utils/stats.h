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
	visited_stats(q);
	//stats on distance comparisons
	distance_stats(q);
}

template<typename T>
auto range_query_stats(parlay::sequence<Tvec_point<T>*> &q){
	auto pred = [&] (Tvec_point<T>* p) {return (p->ngh.size()==0);};
	auto pred1 = [&] (Tvec_point<T>* p) {return !pred(p);};
	auto zero_queries = parlay::filter(q, pred);
	auto nonzero_queries = parlay::filter(q, pred1);
	// std::cout << std::endl;
	// std::cout << "For nonzero entries: " << std::endl;
	parlay::sequence<int> vz = visited_stats(nonzero_queries);
	parlay::sequence<int> dz = distance_stats(nonzero_queries);
	parlay::sequence<int> rz = rounds_stats(nonzero_queries);
	// std::cout << std::endl;
	// std::cout << "For zero entries: " << std::endl;
	parlay::sequence<int> vn = visited_stats(zero_queries);
	parlay::sequence<int> dn = distance_stats(zero_queries);
	parlay::sequence<int> rn = rounds_stats(zero_queries);
	// std::cout << std::endl;
	auto result = {rn, dn, vn, rz, dz, vz};
	return parlay::flatten(result);
}

template<typename T> 
parlay::sequence<int> visited_stats(parlay::sequence<Tvec_point<T>*> &q){
	auto visited_stats = parlay::tabulate(q.size(), [&] (size_t i) {return q[i]->visited;});
	parlay::sort_inplace(visited_stats);
	int avg_visited = (int) parlay::reduce(visited_stats)/((double) q.size());
	size_t tail_index = .99*((float) q.size());
	int tail_visited = visited_stats[tail_index];
	auto result = {avg_visited, tail_visited};
	return result;
	// std::cout << "Average num visited: " << avg_visited << ", 99th percentile num visited: " << tail_visited << std::endl;
}

template<typename T> 
parlay::sequence<int> distance_stats(parlay::sequence<Tvec_point<T>*> &q){
	auto dist_stats = parlay::tabulate(q.size(), [&] (size_t i) {return q[i]->dist_calls;});
	parlay::sort_inplace(dist_stats);
	int avg_dist = (int) parlay::reduce(dist_stats)/((double) q.size());
	size_t tail_index = .99*((float) q.size());
	int tail_dist = dist_stats[tail_index];
	auto result = {avg_dist, tail_dist};
	return result;
	// std::cout << "Average dist cmps: " << avg_dist << ", 99th percentile dist cmps: " << tail_dist << std::endl;
}

template<typename T> 
parlay::sequence<int> rounds_stats(parlay::sequence<Tvec_point<T>*> &q){
	auto exp_stats = parlay::tabulate(q.size(), [&] (size_t i) {return q[i]->rounds;});
	parlay::sort_inplace(exp_stats);
	int avg_exps = (int) parlay::reduce(exp_stats)/((double) q.size());
	size_t tail_index = .99*((float) q.size());
	int tail_exps = exp_stats[tail_index];
	// std::cout << "Average rounds: " << avg_exps << ", 99th percentile rounds: " << tail_exps << 
		// ", max rounds: " << exp_stats[exp_stats.size()-1] << std::endl;
	auto result = {avg_exps, tail_exps, exp_stats[exp_stats.size()-1]};
	return result;
}

void range_gt_stats(parlay::sequence<ivec_point> groundTruth){
  auto sizes = parlay::tabulate(groundTruth.size(), [&] (size_t i) {return groundTruth[i].coordinates.size();});
  parlay::sort_inplace(sizes);
  size_t first_nonzero_index;
  for(size_t i=0; i<sizes.size(); i++){ if(sizes[i] != 0){first_nonzero_index = i; break;}}
  auto nonzero_sizes = (sizes).cut(first_nonzero_index, sizes.size());
  auto sizes_sum = parlay::reduce(nonzero_sizes);
  float avg = static_cast<float>(sizes_sum)/static_cast<float>(nonzero_sizes.size());
  std::cout << "Among nonzero entries, the average number of matches is " << avg << std::endl;
  std::cout << "25th percentile: " << nonzero_sizes[.25*nonzero_sizes.size()] << std::endl;
  std::cout << "75th percentile: " << nonzero_sizes[.75*nonzero_sizes.size()] << std::endl;
  std::cout << "99th percentile: " << nonzero_sizes[.99*nonzero_sizes.size()] << std::endl;
  std::cout << "Max: " << nonzero_sizes[nonzero_sizes.size()-1] << std::endl;
}

