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
// #include "parse_results.h"
#include "types.h"
#include "beamSearch.h"

template<typename T>
range_result checkRecall(
        parlay::sequence<Tvec_point<T>*> &v,
        parlay::sequence<Tvec_point<T>*> &q,
        parlay::sequence<ivec_point> groundTruth,
        int k,
        int beamQ,
        float cut,
        double rad,
        double slack,
        bool random=true,
        int start_point=0) {
  //run twice
  parlay::internal::timer t;
  unsigned d = (v[0]->coordinates).size();
  float query_time;
  if(random){
    rangeSearchRandom(q, v, beamQ, d, rad, k, cut, slack);
    t.next_time();
    rangeSearchRandom(q, v, beamQ, d, rad, k, cut, slack);
    query_time = t.next_time();
  } else{
    rangeSearchAll(q, v, beamQ, d, v[start_point], rad, k, cut, slack);
    t.next_time();
    rangeSearchAll(q, v, beamQ, d, v[start_point], rad, k, cut, slack);
    query_time = t.next_time();
  }
  
  //for range search, disambiguate zero and nonzero queries
  float nonzero_correct = 0.0;
  float zero_correct = 0.0;
  int num_nonzero=0;
  int num_zero=0;

  size_t n = q.size();
  int numCorrect = 0;
  for(int i=0; i<n; i++){
    if(groundTruth[i].coordinates.size() == 0){
      num_zero += 1;
      if(q[i]->ngh.size()==0) {zero_correct += 1.0;}
    }else{
      //since the graph-based nearest neighbor algorithms are exact, no need to check IDs
      num_nonzero += 1;
      int num_real_results = groundTruth[i].coordinates.size();
      int num_correctly_reported = q[i]->ngh.size();
      nonzero_correct += static_cast<float>(num_correctly_reported)/static_cast<float>(num_real_results);
    }
  }
  
  float nonzero_recall = nonzero_correct/static_cast<float>(num_nonzero);
  float zero_recall = zero_correct/static_cast<float>(num_zero);
  float total_recall = (nonzero_correct + zero_correct)/static_cast<float>(num_nonzero + num_zero);

  float QPS = q.size()/query_time;

  auto res = range_query_stats(q);

  range_result R(nonzero_recall, zero_recall, total_recall, res, QPS, k, beamQ, cut, slack);
  return R;
}

template<typename T>
void search_and_parse(parlay::sequence<Tvec_point<T>*> &v, parlay::sequence<Tvec_point<T>*> &q, 
    parlay::sequence<ivec_point> groundTruth, double rad, bool random=true, int start_point=0){
    unsigned d = v[0]->coordinates.size();

    parlay::sequence<range_result> R;
    std::vector<float> slacks = {1.0, 1.5, 2.0, 3.0};
    std::vector<int> beams = {15, 20, 30, 50, 75, 100, 125, 250, 500};
    std::vector<int> allk = {10, 15, 20, 30, 50, 100};
    for(float slack : slacks){
        for(float Q : beams){
            for(float K : allk){
                if(Q>K) R.push_back(checkRecall(v, q, groundTruth, K, Q, 1.14, rad, slack, random, start_point));
            }
        }
    }

    // check "best accuracy"
    R.push_back(checkRecall(v, q, groundTruth, 100, 1000, 10.0, rad, 5.0, random, start_point));

    parlay::sequence<float> buckets = {.4, .5, .6, .7, .8, .9};
    parse_result(R, buckets);
}
