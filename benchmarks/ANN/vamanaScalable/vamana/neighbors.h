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

extern bool report_stats;

template<typename T>
void checkRecall(knn_index<T>& I,
		  parlay::sequence<Tvec_point<T>*> &v,
		  parlay::sequence<Tvec_point<T>*> &q,
		  parlay::sequence<ivec_point> groundTruth,
		  int k,
		  int beamQ,
		  float cut) {
  parlay::internal::timer t;
  int r = 10;
  I.searchNeighbors(q, v, beamQ, k, cut);
  t.next_time();
  I.searchNeighbors(q, v, beamQ, k, cut);
  float query_time = t.next_time();
  float recall = 0.0;
  if (groundTruth.size() > 0) {
    size_t n = q.size();
    int numCorrect = 0;
    for(int i=0; i<n; i++){
      std::set<int> reported_nbhs;
      for(int l=0; l<r; l++)
	reported_nbhs.insert((q[i]->ngh)[l]);
      for(int l=0; l<r; l++)
	if (reported_nbhs.find((groundTruth[i].coordinates)[l])
	    != reported_nbhs.end())
	  numCorrect += 1;
    }
    recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
  }
  std::cout << "k = " << k << ", Q = " << beamQ << ", cut = " << cut
	    << ", throughput = " << (q.size()/query_time) << "/second";
  if (groundTruth.size() > 0)
    std::cout << ", recall = " << recall << std::endl;
  else std::cout << std::endl;
}

template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
	 int beamSize, int beamSizeQ, double alpha, double dummy,
	 parlay::sequence<Tvec_point<T>*> &q,
	 parlay::sequence<ivec_point> groundTruth) {
  parlay::internal::timer t("ANN",report_stats); 
  unsigned d = (v[0]->coordinates).size();
  using findex = knn_index<T>;
  findex I(maxDeg, beamSize, alpha, d);
  parlay::sequence<int> inserts = parlay::tabulate(v.size(), [&] (size_t i){
					    return static_cast<int>(i);});
  I.build_index(v, inserts);
  t.next("Build index");
  std::cout << "num queries = " << q.size() << std::endl;
  std::vector<int> beams = {32, 50, 75, 125, 500};
  std::vector<int> allk = {10, 15, 20, 30, 50, 100};
  std::vector<float> cuts = {1.1, 1.15, 1.2, 1.25};
  for (float cut : cuts)
    for (float Q : beams) 
      checkRecall(I, v, q, groundTruth, 10, Q, cut);

  std::cout << " ... " << std::endl;

  for (float cut : cuts)
    for (int kk : allk)
      checkRecall(I, v, q, groundTruth, kk, 500, cut);

  I.searchNeighbors(q, v, beamSizeQ, k, 1.14);
  if(report_stats){
    graph_stats(v);
    query_stats(q);
    t.next("stats");
  }
}


template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize, double alpha, double dummy) {
  parlay::internal::timer t("ANN",report_stats); 
  { 
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl; 
    using findex = knn_index<T>;
    findex I(maxDeg, beamSize, alpha, d);
    I.build_index(v, parlay::tabulate(v.size(), [&] (size_t i){return static_cast<int>(i);}));
    t.next("Built index");  
    if(report_stats){
      graph_stats(v);
      t.next("stats");
    }
  };
}


    // int parts = 10;
    // size_t m = (size_t) (v.size()/parts);
    // for(int i=0; i<parts; i++){
    //   parlay::sequence<int> delete_list = parlay::tabulate(m, [&] (size_t j){return static_cast<int>(parts*i+j);});
    //   I.lazy_delete(delete_list, v);
    //   I.consolidate_deletes(v);
    //   I.batch_insert(delete_list, v, true);
    // }
