#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/get_time.h" 

#define NOTSELF 1
#include "../include/sfcnn.hpp"
#include "../include/dpoint.hpp"
#include "../KNN/helper.h"


bool report_stats = true;
bool check_correctness = false;
int algorithm_version = 0; 

template<int maxK, class vtx, int Dim>
void ANN_(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN", report_stats);
  uint n = v.size();
  typedef reviver::dpoint<uint, Dim> Point;
  Point *P;
  P = (Point*) parlay::p_malloc(n*sizeof(Point));
  NN_helper<Dim> N;

  N.convert(v, n, P);
  t.next("convert to Kumar's types");

  sfcnn<Point, Dim, uint> NN(P, n);
  t.next("initialize scfnn");

  auto answers = parlay::tabulate(n, [&] (size_t) {
      return std::vector<unsigned long>(k);});
  t.next("initialize answers");

  parlay::parallel_for(0, n, [&] (uint i){
      NN.ksearch(P[i], k, answers[i]);});
  t.next("find all");

  if (check_correctness) {
    parlay::parallel_for(0, n, [&] (uint i){
	if (do_check_correct()) 
	  N.check_correct(i, answers[i][0], P, n);
      });
    t.next("check correctness");	  
  }
}

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
  int d = (v[0]->pt.dimension());
  if (d==2) ANN_<maxK, vtx, 2>(v, k);
  else ANN_<maxK, vtx, 3>(v, k);
}
