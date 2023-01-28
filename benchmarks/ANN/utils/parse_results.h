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

struct range_result{
  double recall;
  double zero_recall;
  double overall_recall;

  int avg_rounds;
  int tail_rounds;
  int max_rounds;

  int avg_cmps;
  int tail_cmps;

  int avg_visited;
  int tail_visited;

  int avg_z_rounds;
  int tail_z_rounds;
  int max_z_rounds;

  int avg_z_cmps;
  int tail_z_cmps;

  int avg_z_visited;
  int tail_z_visited;

  float QPS;

  int k;
  int beamQ;
  float cut;
  double slack;

  range_result(double nzr, double zr, double r, parlay::sequence<int> stats, float qps, int K, int Q, float c, float s) : recall(nzr), zero_recall(zr), 
    overall_recall(r), QPS(qps), k(K), beamQ(Q), cut(c), slack(s) {

    if(stats.size() != 14) abort();

    avg_rounds = stats[0]; tail_rounds = stats[1]; max_rounds = stats[2];
    avg_cmps = stats[3]; tail_cmps = stats[4];
    avg_visited = stats[5]; tail_visited = stats[6];
    avg_z_rounds = stats[7]; tail_z_rounds = stats[8]; max_z_rounds = stats[9];
    avg_z_cmps = stats[10]; tail_z_cmps = stats[11];
    avg_z_visited = stats[12]; tail_z_visited = stats[13];
  }

  void print(){
    std::cout << "k = " << k << ", Q = " << beamQ << ", cut = " << cut << ", slack = " << slack
	    << ", throughput = " << QPS << "/second" << std::endl;
    std::cout << std::endl;
    std::cout << "Nonzero recall: " << recall << std::endl; 
    std::cout << "Zero recall: " << zero_recall << std::endl;
    std::cout << "Combined recall: " << overall_recall << std::endl;
    std::cout << std::endl;
    std::cout << "For nonzero entries: " << std::endl;
    std::cout << "Average rounds: " << avg_rounds << ", 99th percentile rounds: " << tail_rounds << 
		  ", max rounds: " << max_rounds << std::endl;
  	std::cout << "Average dist cmps: " << avg_cmps << ", 99th percentile dist cmps: " << tail_cmps << std::endl;
  	std::cout << "Average num visited: " << avg_visited << ", 99th percentile num visited: " << tail_visited << std::endl;
    std::cout << std::endl; 
    std::cout << "For zero entries: " << std::endl;
    std::cout << "Average rounds: " << avg_z_rounds << ", 99th percentile rounds: " << tail_z_rounds << 
		  ", max rounds: " << max_z_rounds << std::endl;
  	std::cout << "Average dist cmps: " << avg_z_cmps << ", 99th percentile dist cmps: " << tail_z_cmps << std::endl;
  	std::cout << "Average num visited: " << avg_z_visited << ", 99th percentile num visited: " << tail_z_visited << std::endl;
  }
};

struct nn_result{
  double recall;

  int avg_cmps;
  int tail_cmps;

  int avg_visited;
  int tail_visited;

  float QPS;

  int k;
  int beamQ;
  float cut;

  long num_queries;

  nn_result(double r, parlay::sequence<int> stats, float qps, int K, int Q, float c, long q) : recall(r), 
    QPS(qps), k(K), beamQ(Q), cut(c), num_queries(q) {

    if(stats.size() != 4) abort();

    avg_cmps = stats[0]; tail_cmps = stats[1];
    avg_visited = stats[2]; tail_visited = stats[3];
  }

  void print(){
    std::cout << "Over " << num_queries << " queries" << std::endl;
    std::cout << "k = " << k << ", Q = " << beamQ << ", cut = " << cut 
	    << ", throughput = " << QPS << "/second" << std::endl;
    std::cout << "Recall: " << recall << std::endl; 
  	std::cout << "Average dist cmps: " << avg_cmps << ", 99th percentile dist cmps: " << tail_cmps << std::endl;
  	std::cout << "Average num visited: " << avg_visited << ", 99th percentile num visited: " << tail_visited << std::endl;
  }
};

template<typename res>
parlay::sequence<nn_result> parse_result(parlay::sequence<res> results, parlay::sequence<float> buckets){
  parlay::sequence<res> retval;
  for(float b : buckets){
    std::cout << "For recall above: " << b << std::endl;
    auto pred = [&] (res R) {return R.recall >= b;};
    auto candidates = parlay::filter(results, pred);
    if(candidates.size() != 0){
      auto less = [&] (res R, res S) {return R.QPS < S.QPS;};
      res M = *(parlay::max_element(candidates, less));
      M.print();
      retval.push_back(M);
    } else{
      std::cout << "No results found " << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  return retval;
}