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
#include "../../octTree/oct_tree.h"

bool report_stats = false;
#include "../include/parlay_sfcnn_knng.hpp"
#include "../include/dpoint.hpp"
#include "../KNN/helper.h"

bool check_correctness = false;
int algorithm_version = 0; 
int d = 100000; 

template<int maxK, class vtx, int Dim>
void ANN_(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN", report_stats);

  uint n = v.size();
  typedef reviver::dpoint<uint, Dim> Point;       
  Point *P;
  P = (Point*) parlay::p_malloc(n*sizeof(Point)*1);                  
  NN_helper<Dim> N;

  N.convert(v, n, P);
  t.next("convert to Kumar's types");

  sfcnn_knng<Point, Dim, uint> NN(P, n, k); 
  t.next("create nearest neighbor graph");

  if (check_correctness){
    parlay::parallel_for(0, n, [&] (uint i) {
	 if (do_check_correct()){ 
	   N.check_correct(i, NN[i][0], P, n);}});
    t.next("check correctness");
  }
  parlay::p_free(P);
  t.next("delete data");
}

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
  size_t n = v.size(); 
  size_t rounds = n/d; 
  size_t used = 0; 
  for(size_t i = 0; i < rounds; i++){

    parlay::sequence<vtx*> v1;
    v1 = parlay::sequence<vtx*>(d);
    parlay::parallel_for(0, d, [&] (size_t i){
      v1[i] = v[used+i];
    });
    used += d; 
    int d = (v[0]->pt.dimension());
    if (d==2) ANN_<maxK, vtx, 2>(v1, k);
    else ANN_<maxK, vtx, 3>(v1, k);
  }
}
