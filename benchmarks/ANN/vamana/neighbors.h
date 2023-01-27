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
#include "../utils/NSGDist.h"
#include "../utils/types.h"
#include "index.h"
#include "../utils/beamSearch.h"
#include "../utils/indexTools.h"
#include "../utils/stats.h"
#include "../utils/parse_results.h"
#include "../utils/check_nn_recall.h"

extern bool report_stats;

// template<typename T>
// nn_result checkRecall(knn_index<T>& I,
// 		  parlay::sequence<Tvec_point<T>*> &v,
// 		  parlay::sequence<Tvec_point<T>*> &q,
// 		  parlay::sequence<ivec_point> groundTruth,
// 		  int k,
// 		  int beamQ,
// 		  float cut) {
//   parlay::internal::timer t;
//   int r = 10;
//   unsigned d = (v[0]->coordinates).size();
//   I.searchNeighbors(q, v, beamQ, k, cut);
//   t.next_time();
//   I.searchNeighbors(q, v, beamQ, k, cut);
//   float query_time = t.next_time();
//   float recall = 0.0;
//   if (groundTruth.size() > 0) {
//     size_t n = q.size();
//     int numCorrect = 0;
//     for(int i=0; i<n; i++){
//       std::set<int> reported_nbhs;
//       for(int l=0; l<r; l++)
// 	reported_nbhs.insert((q[i]->ngh)[l]);
//       for(int l=0; l<r; l++)
// 	if (reported_nbhs.find((groundTruth[i].coordinates)[l])
// 	    != reported_nbhs.end())
// 	  numCorrect += 1;
//     }
//     recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
//   }
//   float QPS = q.size()/query_time;
//   parlay::sequence<int> stats = query_stats(q);
//   nn_result N(recall, stats, QPS, k, beamQ, cut);
//   return N;
// }

template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
	 int beamSize, int beamSizeQ, double alpha, double dummy,
	 parlay::sequence<Tvec_point<T>*> &q,
	 parlay::sequence<ivec_point> groundTruth, bool graph_built) {
  parlay::internal::timer t("ANN",report_stats);
  unsigned d = (v[0]->coordinates).size();
  using findex = knn_index<T>;
  findex I(maxDeg, beamSize, alpha, d);
  if(graph_built){
    I.find_approx_medoid(v);
  } else{
    parlay::sequence<int> inserts = parlay::tabulate(v.size(), [&] (size_t i){
					    return static_cast<int>(i);});
    I.build_index(v, inserts);
    t.next("Build index");
  }

  int medoid = I.get_medoid();
  search_and_parse(v, q, groundTruth, false, medoid);
  graph_stats(v);
}


template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize, double alpha, double dummy, bool graph_built) {
  parlay::internal::timer t("ANN",report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    using findex = knn_index<T>;
    findex I(maxDeg, beamSize, alpha, d);
    if(graph_built) I.find_approx_medoid(v);
    else{
      parlay::sequence<int> inserts = parlay::tabulate(v.size(), [&] (size_t i){
					    return static_cast<int>(i);});
      I.build_index(v, inserts);
      t.next("Built index");
    }
    if(report_stats){
      graph_stats(v);
      t.next("stats");
    }
  };
}

