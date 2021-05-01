#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include <type_traits> 
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/get_time.h" 

#define NOTSELF 1
#include "../include/parlay_sfcnn.hpp"
#include "../include/dpoint.hpp"
#include "../KNN/helper.h"

int d = 20000000; 
bool report_stats = false;
bool check_correctness = false;
int algorithm_version = 0;

template<int maxK, class vtx, int Dim>
void ANN_(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN", report_stats);{
    size_t rounds = v.size()/d; 
    std::cout << rounds << " rounds" << std::endl; 
    size_t used = 0; 
    for(size_t i = 0; i < rounds; i++){

      parlay::sequence<vtx*> v1;
      v1 = parlay::sequence<vtx*>(d);
      parlay::parallel_for(0, d, [&] (size_t i){
        v1[i] = v[used+i];
      });
      used += d; 


      uint n = v1.size();
      uint h = n/2; 
      typedef reviver::dpoint<uint, Dim> Point;
      Point *P;
      P = (Point*) parlay::p_malloc(h*sizeof(Point));
      Point *Q;
      Q = (Point*) parlay::p_malloc(h*sizeof(Point));
      NN_helper<Dim, vtx> N;

      N.convert(v1, n, P, Q);
      t.next("convert to Kumar's types");     

      sfcnn<Point, Dim, uint> NN(P, h);
      t.next("initialize scfnn");

      auto answers = parlay::tabulate(h, [&] (size_t) {
          return parlay::sequence<unsigned long>(k);});
      t.next("initialize answers");

      parlay::parallel_for(0, h, [&] (uint i){
          NN.ksearch(Q[i], k, answers[i]);});
      t.next("find all");

      // if (check_correctness) {
      //   parlay::parallel_for(0, n, [&] (uint i){
      // if (do_check_correct()) 
      //   N.check_correct(i, answers[i][0], P, n);
      //     });
      //   t.next("check correctness");    
      // }
  }
  }
}

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
  int d = (v[0]->pt.dimension());
  if (d==2) ANN_<maxK, vtx, 2>(v, k);
  else ANN_<maxK, vtx, 3>(v, k);
}

