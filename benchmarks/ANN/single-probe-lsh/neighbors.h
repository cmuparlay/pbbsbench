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
#include "../utils/beamSearch.h"
#include "../utils/indexTools.h"
#include "../utils/stats.h"
#include "../utils/parse_results.h"

#include "lsh.h"
#include "../utils/csvfile.h"

extern bool report_stats;


// Going to simply duplicate code from check_nn_recall.h for now...
namespace recall {

template <typename T, typename I>
void searchAll(I& index, parlay::sequence<Tvec_point<T>*>& q, int k, unsigned dim,
  parlay::sequence<uint32_t>& query_result_size,
  parlay::sequence<uint32_t>& query_result_ids,
  parlay::sequence<float>&    query_result_dists) {

  grann::Parameters search_params;

  std::cout << "Running searches" << std::endl;
	parlay::parallel_for(0, q.size(), [&] (size_t i) {
    uint32_t res_count = index.search(q[i]->coordinates.begin(), dim, k, search_params,
        query_result_ids.begin() + i * k,
        query_result_dists.begin() + i * k);
    query_result_size[i] = res_count;
  });
  std::cout << "Finished searches" << std::endl;

}

template<typename T, typename I>
nn_result checkRecall(
        I& index,
        parlay::sequence<Tvec_point<T>*> &q,
        parlay::sequence<ivec_point> groundTruth,
        int k,
        unsigned d) {

  parlay::sequence<uint32_t>    query_result_size(q.size());
  parlay::sequence<uint32_t> query_result_ids(q.size() * k);
  parlay::sequence<float>    query_result_dists(q.size() * k);

  parlay::internal::timer t;
  int r = 10;
  float query_time;

  searchAll(index, q, k, d, query_result_size, query_result_ids, query_result_dists);
  t.next_time();
  searchAll(index, q, k, d, query_result_size, query_result_ids, query_result_dists);
  query_time = t.next_time();

  for (size_t i=0; i<q.size(); ++i) {
    q[i]->ngh.resize(query_result_size[i]);
    for (size_t j=0; j<query_result_size[i]; ++j) {
      q[i]->ngh[j] = query_result_ids[i*k + j];
    }
  }

  float recall = 0.0;
  bool dists_present = (groundTruth[0].distances.size() != 0);
  if (groundTruth.size() > 0 && !dists_present) {
    size_t n = q.size();
    int numCorrect = 0;
    for(int i=0; i<n; i++){
      std::set<int> reported_nbhs;
      for(int l=0; l<r; l++) reported_nbhs.insert((q[i]->ngh)[l]);
      for(int l=0; l<r; l++){
	      if (reported_nbhs.find((groundTruth[i].coordinates)[l]) != reported_nbhs.end()){
          numCorrect += 1;
      }
      }
    }
    recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
  } else if(groundTruth.size() > 0 && dists_present) {
    size_t n = q.size();
    int numCorrect = 0;
    for(int i=0; i<n; i++){
      parlay::sequence<int> results_with_ties;
      for(int l=0; l<r; l++) results_with_ties.push_back(groundTruth[i].coordinates[l]);
      float last_dist = groundTruth[i].distances[r-1];
      for(int l=r; l<groundTruth[i].coordinates.size(); l++){
        if(groundTruth[i].distances[l] == last_dist){
          results_with_ties.push_back(groundTruth[i].coordinates[l]);
        }
      }
      std::set<int> reported_nbhs;
      for(int l=0; l<r; l++) reported_nbhs.insert((q[i]->ngh)[l]);

      for(int l=0; l<results_with_ties.size(); l++){
	      if (reported_nbhs.find(results_with_ties[l]) != reported_nbhs.end()){
          numCorrect += 1;
      }
      }
    }
    recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
  }
  float QPS = q.size()/query_time;
  auto stats = query_stats(q);
  std::cout << "recall = " << recall << " QPS = " << QPS << " k = " << k << std::endl;
  // Q, c are just some dummy values.
  nn_result N(recall, stats, QPS, k, 500, 500, q.size());
  return N;
}

void write_to_csv(std::string csv_filename, parlay::sequence<float> buckets,
  parlay::sequence<nn_result> results){
  csvfile csv(csv_filename);
  //csv << "GRAPH" << "Parameters" << "Size" << "Build time" << "Avg degree" << "Max degree" << endrow;
  //csv << G.name << G.params << G.size << G.time << G.avg_deg << G.max_deg << endrow;
  csv << endrow;
  csv << "Num queries" << "Target recall" << "Actual recall" << "QPS" << "Average Cmps" <<
    "Tail Cmps" << "Average Visited" << "Tail Visited" << "k" << "Q" << "cut" << endrow;
  for(int i=0; i<results.size(); i++){
    nn_result N = results[i];
    csv << N.num_queries << buckets[i] << N.recall << N.QPS << N.avg_cmps << N.tail_cmps <<
      N.avg_visited << N.tail_visited << N.k << N.beamQ << N.cut << endrow;
  }
  csv << endrow;
  csv << endrow;
}

template<typename T, typename I>
void search_and_parse(I& index, parlay::sequence<Tvec_point<T>*> &v, parlay::sequence<Tvec_point<T>*> &q,
    parlay::sequence<ivec_point> groundTruth, char* res_file) {
    unsigned d = v[0]->coordinates.size();

    parlay::sequence<nn_result> results;
    std::vector<int> beams = {15, 20, 30, 50, 75, 100, 125, 250, 500};
    std::vector<int> allk = {10, 15, 20, 30, 50, 100};
    std::vector<float> cuts = {1.1, 1.125, 1.15, 1.175, 1.2, 1.25};

    for (int kk : allk)
      results.push_back(checkRecall(index, q, groundTruth, kk, d));

    parlay::sequence<float> buckets = {.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .73, .75, .77, .8, .83, .85, .87, .9, .93, .95, .97, .99, .995, .999};
    auto [res, ret_buckets] = parse_result(results, buckets);
    if(res_file != NULL) write_to_csv(std::string(res_file), ret_buckets, res);
}

}  // namespace recall


template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
	 int beamSize, int beamSizeQ, double alpha, double dummy,
	 parlay::sequence<Tvec_point<T>*> &q,
	 parlay::sequence<ivec_point> groundTruth, char* res_file, bool graph_built) {
  std::cout << "Here" << std::endl;

  grann::LSHIndex<T> I(v, grann::L2);

  grann::Parameters params;
  // Hackily overload graph-based parameters:
  params.Set<uint32_t>("num_tables", maxDeg);
  params.Set<uint32_t>("table_size", beamSize);

  parlay::internal::timer t("ANN",report_stats);
  I.build(params);
  double idx_time = t.next_time();

  // std::string name = "Single-Probe-LSH";
  // std::string params = "NumTables = " + std::to_string(maxDeg) + ", TableSize = " + std::to_string(beamSize);

  parlay::sequence<uint32_t> query_result_ids(q.size() * k);
  parlay::sequence<float>    query_result_dists(q.size() * k);

  uint32_t dim = q[0]->coordinates.size();

  recall::search_and_parse(I, v, q, groundTruth, res_file);
}


template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize, double alpha, double dummy, bool graph_built) {
}

