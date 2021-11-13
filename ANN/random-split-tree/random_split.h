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
#include <random>
#include <set>

template<class fvec_point>
struct random_split_tree {

	using slice_fvec = decltype(make_slice(parlay::sequence<fvec_point*>()));

	random_split_tree() {}

  std::pair<size_t, size_t> select_two_random(
      parlay::sequence<fvec_point*>& points,
      parlay::sequence<size_t>& active_indices,
      parlay::random& rnd) {
    size_t first_index = rnd.ith_rand(0) % active_indices.size(); 
    size_t second_index_unshifted = rnd.ith_rand(1) % (active_indices.size()-1);
    size_t second_index = (second_index_unshifted < first_index) ?
      second_index_unshifted : (second_index_unshifted + 1);
//    std::cout << "Active indices size = " << active_indices.size() <<
//      " fst = " << first_index << " snd = " << second_index <<
//      std::endl;

    return {active_indices[first_index], active_indices[second_index]};
  }

  // Returns the depth of the tree.
  size_t random_topdown_cluster(parlay::sequence<fvec_point*>& points,
      parlay::sequence<size_t>& active_indices, parlay::random& rnd) {

    if (active_indices.size() < 20) {
//      std::cout << "Done, stopping at size = " << active_indices.size() << std::endl;
      return 1;
    }

//    std::cout << "Cluster size = " << active_indices.size() << std::endl;
    // select two random points
    auto [u, v] = select_two_random(points, active_indices, rnd);
    fvec_point* first = points[u];
    fvec_point* second = points[v];
//    std::cout << "u = " << u << " v = " << v << std::endl;

    // Split points based on which of the two points are closer.
    auto closer_first = parlay::filter(parlay::make_slice(active_indices), [&] (size_t ind) {
      fvec_point* p = points[ind];
      auto dist_first = distance(p, first);
      auto dist_second = distance(p, second);
      if (ind != u && ind != v) {
        return dist_first <= dist_second;
      }
      return false;
    });

    auto closer_second = parlay::filter(parlay::make_slice(active_indices), [&] (size_t ind) {
      fvec_point* p = points[ind];
      auto dist_first = distance(p, first);
      auto dist_second = distance(p, second);
      if (ind != u && ind != v) {
        return dist_second < dist_first;
      }
      return false;
    });

    auto left_rnd = rnd.fork(0);
    auto right_rnd = rnd.fork(1);
    size_t left_depth = 0;
    size_t right_depth = 0;
    parlay::par_do([&] () { left_depth = random_topdown_cluster(points, closer_second, left_rnd); },
      [&] () { right_depth = random_topdown_cluster(points, closer_first, right_rnd); });

    return 1 + std::max(left_depth, right_depth);
	}

	void build_index(parlay::sequence<fvec_point*> v) {
    std::random_device rd;    
  	std::mt19937 rng(rd());   
  	std::uniform_int_distribution<int> uni(0,v.size()); 


    parlay::random rnd(uni(rng));
    auto active_indices = parlay::tabulate(v.size(), [&] (size_t i) { return i; });
		size_t max_depth = random_topdown_cluster(v, active_indices, rnd);
    std::cout << "Max Depth = " << max_depth << std::endl;
	}

};
