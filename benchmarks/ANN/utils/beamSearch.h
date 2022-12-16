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

#ifndef BEAMSEARCH
#define BEAMSEARCH

#include <algorithm>
#include <set>
#include <unordered_set>
#include "common/geometry.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "types.h"
#include "indexTools.h"

extern bool report_stats;

using pid = std::pair<int, float>;

// returns true if F \setminus V = emptyset
bool intersect_nonempty(parlay::sequence<pid>& V, parlay::sequence<pid>& F) {
  for (int i = 0; i < F.size(); i++) {
    auto pred = [&](pid a) { return F[i].first == a.first; };
    if (parlay::find_if(V, pred) == V.end()) return true;
  }
  return false;
}

// will only be used when there is an element in F that is not in V
// hence the ``return 0" line will never be called
pid id_next(parlay::sequence<pid>& V, parlay::sequence<pid>& F) {
  for (int i = 0; i < F.size(); i++) {
    auto pred = [&](pid a) { return F[i].first == a.first; };
    if (parlay::find_if(V, pred) == V.end()) return F[i];
  }
  return std::make_pair(0, 0);
}

// for debugging
void print_seq(parlay::sequence<int> seq) {
  int fsize = seq.size();
  std::cout << "[";
  for (int i = 0; i < fsize; i++) {
    std::cout << seq[i] << ", ";
  }
  std::cout << "]" << std::endl;
}

// updated version by Guy
template <typename T>
std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
    Tvec_point<T>* p, parlay::sequence<Tvec_point<T>*>& v,
    Tvec_point<T>* starting_point, int beamSize, unsigned d, int k=0, float cut=1.14) {
  
  parlay::sequence<Tvec_point<T>*> start_points;
  start_points.push_back(starting_point);
  return beam_search(p, v, start_points, beamSize, d, k, cut);

}

// updated version by Guy
template <typename T>
std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
    Tvec_point<T>* p, parlay::sequence<Tvec_point<T>*>& v,
    parlay::sequence<Tvec_point<T>*> starting_points, int beamSize, unsigned d, int k=0, float cut=1.14) {
  // initialize data structures
  auto vvc = v[0]->coordinates.begin();
  long stride = v[1]->coordinates.begin() - v[0]->coordinates.begin();
  std::vector<pid> visited;
  auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first); };
  auto make_pid = [&] (int q) {
		    return std::pair{q, distance(vvc + q*stride, p->coordinates.begin(), d)};};
  int bits = std::ceil(std::log2(beamSize*beamSize))-2;
  parlay::sequence<int> hash_table(1 << bits, -1);

  auto pre_frontier = parlay::tabulate(starting_points.size(), [&] (size_t i) {
    return make_pid(starting_points[i]->id);
  });

  auto frontier = parlay::sort(pre_frontier, less);

  std::vector<pid> unvisited_frontier(beamSize);
  parlay::sequence<pid> new_frontier(beamSize + v[0]->out_nbh.size());
  unvisited_frontier[0] = frontier[0];
  int remain = 1;

  // terminate beam search when the entire frontier has been visited
  while (remain > 0) {
    // the next node to visit is the unvisited frontier node that is closest to p
    pid currentPid = unvisited_frontier[0];
    Tvec_point<T>* current = v[currentPid.first];
    auto nbh = current->out_nbh.cut(0, size_of(current->out_nbh));
    auto candidates = parlay::filter(nbh,  [&] (int a) {
	     int loc = parlay::hash64_2(a) & ((1 << bits) - 1);
	     if (a == p->id || hash_table[loc] == a) return false;
	     hash_table[loc] = a;
	     return true;});
    auto pairCandidates = parlay::map(candidates, [&] (long c) {return make_pid(c);}, 1000);
    auto sortedCandidates = parlay::sort(pairCandidates, less);
    auto f_iter = std::set_union(frontier.begin(), frontier.end(),
				 sortedCandidates.begin(), sortedCandidates.end(),
				 new_frontier.begin(), less);
    size_t f_size = std::min<size_t>(beamSize, f_iter - new_frontier.begin());
    if (k > 0 && f_size > k) 
      f_size = (std::upper_bound(new_frontier.begin(), new_frontier.begin() + f_size,
				std::pair{0, cut * new_frontier[k].second}, less)
		- new_frontier.begin());
    frontier = parlay::tabulate(f_size, [&] (long i) {return new_frontier[i];});
    visited.insert(std::upper_bound(visited.begin(), visited.end(), currentPid, less), currentPid);
    auto uf_iter = std::set_difference(frontier.begin(), frontier.end(),
				 visited.begin(), visited.end(),
				 unvisited_frontier.begin(), less);
    remain = uf_iter - unvisited_frontier.begin();
  }
  return std::make_pair(frontier, parlay::to_sequence(visited));
}


// searches every element in q starting from a randomly selected point
template <typename T>
int beamSearchRandom(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      unsigned d, double cut = 1.14) {
  if ((k + 1) > beamSizeQ) {
    std::cout << "Error: beam search parameter Q = " << beamSizeQ
              << " same size or smaller than k = " << k << std::endl;
    abort();
  }
  if(report_stats) efanna2e::distance_calls.reset();
  // use a random shuffle to generate random starting points for each query
  size_t n = v.size();
  auto indices =
      parlay::random_permutation<int>(static_cast<int>(n), time(NULL));
  // for(size_t i=0; i<q.size(); i++){
  parlay::parallel_for(0, q.size(), [&](size_t i) {
    parlay::sequence<int> neighbors = parlay::sequence<int>(k);
    size_t index = indices[i];
    Tvec_point<T>* start = v[index];
    parlay::sequence<pid> beamElts;
    parlay::sequence<pid> visitedElts;
    std::pair<parlay::sequence<pid>, parlay::sequence<pid>> pairElts;
    pairElts = beam_search(q[i], v, start, beamSizeQ, d, k, cut);
    beamElts = pairElts.first;
    visitedElts = pairElts.second;
    // the first element of the frontier may be the point itself
    // if this occurs, do not report it as a neighbor
    if (beamElts[0].first == i) {
      for (int j = 0; j < k; j++) {
        neighbors[j] = beamElts[j + 1].first;
      }
    } else {
      for (int j = 0; j < k; j++) {
        neighbors[j] = beamElts[j].first;
      }
    }
    q[i]->ngh = neighbors;
    if (report_stats) q[i]->cnt = visitedElts.size();
  });
  if(report_stats) return efanna2e::distance_calls.get_value();
  else return 0;
  // }
}

template <typename T>
int searchAll(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      unsigned d, Tvec_point<T>* starting_point, float cut) {
    
    parlay::sequence<Tvec_point<T>*> start_points;
    start_points.push_back(starting_point);
    return searchAll(q, v, beamSizeQ, k, d, start_points, cut);
}

template <typename T>
int searchAll(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      unsigned d, parlay::sequence<Tvec_point<T>*> starting_points, float cut) {
  if ((k + 1) > beamSizeQ) {
    std::cout << "Error: beam search parameter Q = " << beamSizeQ
              << " same size or smaller than k = " << k << std::endl;
    abort();
  }
  if(report_stats) efanna2e::distance_calls.reset();
  parlay::parallel_for(0, q.size(), [&](size_t i) {
    parlay::sequence<int> neighbors = parlay::sequence<int>(k);
    parlay::sequence<pid> beamElts;
    parlay::sequence<pid> visitedElts;
    std::pair<parlay::sequence<pid>, parlay::sequence<pid>> pairElts;
    pairElts = beam_search(q[i], v, starting_points, beamSizeQ, d, 0, cut);
    beamElts = pairElts.first;
    visitedElts = pairElts.second;
    // the first element of the frontier may be the point itself
    // if this occurs, do not report it as a neighbor
    if (beamElts[0].first == i) {
      for (int j = 0; j < k; j++) {
        neighbors[j] = beamElts[j + 1].first;
      }
    } else {
      for (int j = 0; j < k; j++) {
        neighbors[j] = beamElts[j].first;
      }
    }
    q[i]->ngh = neighbors;
    if (report_stats) q[i]->cnt = visitedElts.size();
  });
  if(report_stats) return efanna2e::distance_calls.get_value();
  return 0;
}

template<typename T>
void rangeSearchAll(parlay::sequence<Tvec_point<T>*> q, parlay::sequence<Tvec_point<T>*>& v, 
  int beamSize, unsigned d, Tvec_point<T>* start_point, double r, double slack = 1.0){

    parlay::parallel_for(0, q.size(), [&] (size_t i){
      auto in_range = range_search(q[i], v, beamSize, d, start_point, r, slack);
      parlay::sequence<int> nbh;
      for(auto j : in_range) nbh.push_back(j);
      q[i]->ngh = nbh;
    });
}

template<typename T>
std::set<int> rangeSearch(Tvec_point<T>* q, parlay::sequence<Tvec_point<T>*>& v, 
  int beamSize, unsigned d, Tvec_point<T>* start_point, double r, double slack = 1.0){
  
  double max_rad = 0;

  std::set<int> nbh;

  int beam = beamSize;

  while(max_rad <= slack*r){
    auto [visited, neighbors] = beam_search(q, v, start_point, beam, d);
    double new_rad = neighbors[neighbors.size()-1].second;
    for(auto p : neighbors){
      if(p.second <= r) nbh.insert(p.first);
    }
    beam *= 2;
    max_rad = new_rad;
  }

  return nbh;

}
                      

#endif

