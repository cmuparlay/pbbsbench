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
#include "vamana/index.h"
#include "../utils/beamSearch.h"
#include "../utils/indexTools.h"
#include "../utils/stats.h"
#include "../utils/parse_results.h"
#include "../utils/check_range_recall.h"

extern bool report_stats;


template<typename T>
void RNG(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
	 int beamSize, int beamSizeQ, double alpha, double dummy, double rad,
	 parlay::sequence<Tvec_point<T>*> &q,
	 parlay::sequence<ivec_point> groundTruth, bool graph_built) {
  parlay::internal::timer t("ANN",report_stats);
  // gt_stats(groundTruth);
  unsigned d = (v[0]->coordinates).size();
  using findex = knn_index<T>;
  findex I(maxDeg, beamSize, alpha, d);
  if(graph_built){
    I.find_approx_medoid(v);
    t.next("Find medoid");
  } else{
    parlay::sequence<int> inserts = parlay::tabulate(v.size(), [&] (size_t i){
					    return static_cast<int>(i);});
    I.build_index(v, inserts);
    t.next("Build index");
  }
  int medoid = I.get_medoid();
  search_and_parse(v, q, groundTruth, rad, false, medoid);
  t.next("Searching");
  graph_stats(v);
  t.next("stats");

}
