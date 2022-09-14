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
#include "clustertree.h"
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

// takes in two sorted sequences and returns a sorted union
template <typename Seq>
parlay::sequence<pid> seq_union(Seq& P,
                                parlay::sequence<pid>& Q) {
  auto less = [&](pid a, pid b) { return a.second < b.second; };
  pid* first1 = P.begin();
  pid* last1 = P.end();
  pid* first2 = Q.begin();
  pid* last2 = Q.end();
  parlay::sequence<pid> result = parlay::sequence<pid>();
  result.reserve(P.size()+Q.size());
  while (true) {
    if (first1 == last1) {
      while (first2 != last2) {
        result.push_back(*first2);
        ++first2;
      }
      return result;
    } else if (first2 == last2) {
      while (first1 != last1) {
        result.push_back(*first1);
        ++first1;
      }
      return result;
    }
    if (less(*first1, *first2)) {
      result.push_back(*first1);
      ++first1;
    } else if (less(*first2, *first1)) {
      result.push_back(*first2);
      ++first2;
    } else {
      if (first1->first == first2->first) {
        result.push_back(*first1);
        ++first1;
        ++first2;
      } else {
        result.push_back(*first1);
        result.push_back(*first2);
        ++first1;
        ++first2;
      }
    }
  }
  return result;
}

// updated version by Guy
template <typename T>
std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
    Tvec_point<T>* p, parlay::sequence<Tvec_point<T>*>& v,
    Tvec_point<T>* medoid, int beamSize, unsigned d, int k=0) {
  // initialize data structures
  auto vvc = v[0]->coordinates.begin();
  long stride = v[1]->coordinates.begin() - v[0]->coordinates.begin();
  std::vector<pid> visited;
  parlay::sequence<pid> frontier;
  auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first); };
  auto make_pid = [&] (int q) {
		    return std::pair{q, distance(vvc + q*stride, p->coordinates.begin(), d)};};
  int bits = std::ceil(std::log2(beamSize*beamSize))-2;
  parlay::sequence<int> hash_table(1 << bits, -1);

  // the frontier starts with the medoid
  frontier.push_back(make_pid(medoid->id));

  std::vector<pid> unvisited_frontier(beamSize);
  parlay::sequence<pid> new_frontier(2*beamSize);
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
				std::pair{0, 1.14 * new_frontier[k].second}, less)
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

template <typename T>
std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search_prev(
    Tvec_point<T>* p, parlay::sequence<Tvec_point<T>*>& v,
    Tvec_point<T>* medoid, int beamSize, unsigned d) {
  // initialize data structures
  parlay::sequence<pid> visited = parlay::sequence<pid>();
  visited.reserve(5*beamSize);
  parlay::sequence<pid> frontier = parlay::sequence<pid>();
  frontier.reserve(2*beamSize);
  auto less = [&](pid a, pid b) { return a.second < b.second; };  
  // the frontier starts with the medoid
  frontier.push_back(std::make_pair(
      medoid->id,
      distance(medoid->coordinates.begin(), p->coordinates.begin(), d)));
  // terminate beam search when the entire frontier has been visited
  while (intersect_nonempty(visited, frontier)) {
    // the next node to visit is the unvisited frontier node that is closest to
    // p
    pid currentPid = (id_next(visited, frontier));
    Tvec_point<T>* current = v[currentPid.first];
    auto g = [&](int a) {
      auto pred = [&](pid p) { return p.first == a; };
      if (parlay::find_if(visited, pred) == visited.end()) return true;
      return false;
    };
    parlay::sequence<int> candidatesTemp = parlay::tabulate(size_of(current->out_nbh), [&] (size_t i) {
      return current->out_nbh[i];});
    auto candidates = parlay::filter(candidatesTemp, g);
    parlay::sequence<pid> pairCandidates =
        parlay::sequence<pid>(candidates.size());
    for (int i = 0; i < candidates.size(); i++) {
      pairCandidates[i] = std::make_pair(
          candidates[i], distance(v[candidates[i]]->coordinates.begin(),
                                  p->coordinates.begin(), d));
    }
    auto sortedCandidates = parlay::sort(pairCandidates, less);
    frontier = seq_union(frontier, sortedCandidates);
    if (frontier.size() > beamSize)
      frontier.erase(frontier.begin() + beamSize, frontier.end());
    // add the node to the visited list
    visited.push_back(currentPid);

  }
  return std::make_pair(frontier, visited);
}

// searches every element in q starting from a randomly selected point
template <typename T>
void beamSearchRandom(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      unsigned d) {
  if ((k + 1) > beamSizeQ) {
    std::cout << "Error: beam search parameter Q = " << beamSizeQ
              << " same size or smaller than k = " << k << std::endl;
    abort();
  }
  // use a random shuffle to generate random starting points for each query
  size_t n = v.size();
  auto indices =
      parlay::random_permutation<int>(static_cast<int>(n), time(NULL));
  parlay::parallel_for(0, q.size(), [&](size_t i) {
    parlay::sequence<int> neighbors = parlay::sequence<int>(k);
    Tvec_point<T>* start = v[indices[i]];
    parlay::sequence<pid> beamElts;
    parlay::sequence<pid> visitedElts;
    std::pair<parlay::sequence<pid>, parlay::sequence<pid>> pairElts;
    pairElts = beam_search(q[i], v, start, beamSizeQ, d);
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
}

template <typename T>
void searchFromSingle(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      unsigned d, Tvec_point<T>* medoid) {
  if ((k + 1) > beamSizeQ) {
    std::cout << "Error: beam search parameter Q = " << beamSizeQ
              << " same size or smaller than k = " << k << std::endl;
    abort();
  }
  parlay::parallel_for(0, q.size(), [&](size_t i) {
    parlay::sequence<int> neighbors = parlay::sequence<int>(k);
    parlay::sequence<pid> beamElts;
    parlay::sequence<pid> visitedElts;
    std::pair<parlay::sequence<pid>, parlay::sequence<pid>> pairElts;
    pairElts = beam_search(q[i], v, medoid, beamSizeQ, d, k);
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
}

//UNDER CONSTRUCTION
template<typename T>
void warmStartSearch(parlay::sequence<Tvec_point<T>*>& q,
                      parlay::sequence<Tvec_point<T>*>& v, int beamSizeQ, int k,
                      int cluster_size, int num_clusterings, int cluster_size_inner, unsigned d){
  if(cluster_size_inner < 1){
    std::cout << "Error: cluster size must be larger than 1" << std::endl; abort(); 
  }
  if(num_clusterings*cluster_size > v.size()){
    std::cout << "Error: cannot have more points in cluster tree than contained in dataset" << std::endl;
    abort();
  }
  using ctree = cluster_tree<T>;
  using node = typename ctree::node;
  parlay::sequence<node*> clusters(num_clusterings);
  auto indices = parlay::random_permutation<int>(static_cast<int>(v.size()), time(NULL));
  parlay::parallel_for(0, num_clusterings, [&] (size_t i){
    auto tree_members = parlay::tabulate(cluster_size, [&] (size_t j) {return v[indices[cluster_size*i+j]];});
    clusters[i] = ctree::random_clustering_wrapper(tree_members, cluster_size_inner, d);
  });
  parlay::parallel_for(0, q.size(), [&] (size_t i){
    auto candidate_starts = parlay::tabulate(num_clusterings, [&] (size_t j) {
      return ctree::search(clusters[j], q[i], v, d);
    });
    auto comp = [&] (std::pair<size_t, float> a, std::pair<size_t, float> b) {return a.second < b.second;};
    size_t start =  (*(parlay::min_element(candidate_starts, comp))).first;
    for(const std::pair<size_t, float> i : candidate_starts) std::cout << i.first << std::endl; 
    parlay::sequence<int> neighbors = parlay::sequence<int>(k);
    parlay::sequence<pid> beamElts;
    parlay::sequence<pid> visitedElts;
    std::pair<parlay::sequence<pid>, parlay::sequence<pid>> pairElts;
    pairElts = beam_search(q[i], v, v[start], beamSizeQ, d);
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
}

#endif

// // Uses a std::set for the frontier and beam, and unordered_set for
// // visited elements.
// // NOT CURRENTLY USED ANYWHERE
// template <typename T>
// parlay::sequence<pid> beam_search_2(Tvec_point<T>* p,
//                                     parlay::sequence<Tvec_point<T>*>& v,
//                                     Tvec_point<T>* medoid, int beamSize,
//                                     unsigned d) {
//   // initialize data structures
//   std::unordered_set<int> visited;
//   auto q_comp = [] (const pid& l, const pid& r) {
//     if (l.second < r.second) return true;
//     if (l.second == r.second) return l.first < r.first;
//     return false;
//   };
//   std::set<pid, decltype(q_comp)> frontier(q_comp);
//   std::set<pid, decltype(q_comp)> beam(q_comp);

//   auto ret = parlay::sequence<pid>();
//   ret.reserve(beamSize);

//   // the frontier starts with the medoid
//   auto medoid_pair = std::make_pair(
//       medoid->id,
//       distance(medoid->coordinates.begin(), p->coordinates.begin(), d));
//   frontier.insert(medoid_pair);
//   beam.insert(medoid_pair);

//   // terminate beam search when the entire frontier has been visited
//   while (!frontier.empty()) {
//     // the next node to visit is the unvisited frontier node that is closest to
//     // p
//     pid currentPid = *(frontier.begin());
//     frontier.erase(frontier.begin());

//     ret.push_back(currentPid);
//     visited.insert(currentPid.first);

//     Tvec_point<T>* current = v[currentPid.first];
//     auto g = [&] (int ngh) {
//       if (visited.find(ngh) != visited.end()) return false;
//       return true;
//     };

//     // for (const int& ngh : current->out_nbh) {
//     //   if (visited.find(ngh) == visited.end()) {
//     //     float dist = distance(v[ngh]->coordinates.begin(), p->coordinates.begin(), d);
//     //     // insert it into the beam
//     //     beam.insert(std::make_pair(ngh, dist));
//     //     frontier.insert(std::make_pair(ngh, dist));
//     //     visited.insert(ngh);
//     //   }
//     // }

//     for(int i=0; i<size_of(current->out_nbh); i++){
//       int ngh = current->out_nbh[i];
//       if (visited.find(ngh) == visited.end()) {
//         float dist = distance(v[ngh]->coordinates.begin(), p->coordinates.begin(), d);
//         // insert it into the beam
//         beam.insert(std::make_pair(ngh, dist));
//         frontier.insert(std::make_pair(ngh, dist));
//         visited.insert(ngh);
//       }
//     }

//     while (beam.size() > beamSize) {
//       auto it = std::prev(beam.end());
//       pid point = *it;
//       beam.erase(it);
//       auto f_it = frontier.find(point);
//       if (f_it != frontier.end()) {
//         frontier.erase(f_it);
//       }
//     }
//   }
//   return ret;
// }
