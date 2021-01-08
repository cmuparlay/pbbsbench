#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include <type_traits> 
#include "parlay/parallel.h"
#include "parlay/primitives.h"
// #include "parlay/alloc.h"
#include "common/geometry.h"
#include "common/get_time.h" 
#include "../../octTree/oct_tree.h"
bool parallel = true;
#include "../include/sfcnn.hpp"
#include "../include/dpoint.hpp"
#define sq(x) (((double) (x))* ((double) (x)))


int algorithm_version = 0; 
int key_bits = 64;
double eps = 0;
bool report_stats = true;
bool check_correctness = false;

uint abs(uint x, uint y){
	if (x > y){
		return (x-y);
	} else{
		return (y-x);
	}
}

template<int d>
struct NN_helper{

	typedef reviver::dpoint<uint, d> Point;

	double squared_distance(Point p, Point q){ 
		double dist = 0; 
		for(int j = 0; j<d; j++){
			dist += sq(abs(q[j], p[j]));
		}
		return dist; 
	}

	bool are_equal(Point p, Point q){
		for (int j = 0; j<d; j++){
			if(p[j] != q[j]){
				return false;
			}
		}
		return true;
	}

	void print_point(Point p){
		std::cout << "Point: ";
		for(int j=0; j<d; j++){
			std::cout << p[j] << ", ";
		}
		std::cout << "\n";
	}

	template<class vtx>
	void convert(parlay::sequence<vtx*> &v, uint n, Point P[]){
		//prelims for rounding each point to an integer: 
		// 1) find the smallest point in each dimension
		// 2) find the largest gap between min and max over all dimensions
		using o_tree = oct_tree<vtx>;
		using box = typename o_tree::box;
		using point = typename vtx::pointT;
		box b = o_tree::get_box(v);
		double Delta = 0;
		for (int i = 0; i < d; i++) 
		  Delta = std::max(Delta, b.second[i] - b.first[i]); 
		point min_point = b.first;
		// round each point to an integer with key_bits bit
		int bits = 31; //key_bits/d;
		int maxval = (((size_t) 1) << bits) - 1;
	    parlay::parallel_for(0, n, [&] (uint i){
	    	Point p; 
	    	P[i] = p;
	    	for (int j = 0; j < d; j++){
	      		uint coord = (uint) floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta); 
	      		P[i][j] = coord;
	      	}
	    });
	}

	void check_correct(unsigned long idx, unsigned long idx_other, Point P[], uint n){
	  if (idx == idx_other) {
	    std::cout << "Error: nearest neighbors is self at index " << idx << endl;
	    abort();
	  }
	  Point nearest;
	  double nearest_dist = DBL_MAX;
	  Point q = P[idx];
	  double reported_distance = squared_distance(P[idx_other], q);
	  for(uint i = 0; i < n; i++){
	    if(i != idx){ //make sure we don't report the query point as its own nearest neighbor
	      double dist = squared_distance(q, P[i]);
	      if(dist < nearest_dist){
		nearest = P[i];
		nearest_dist = dist; 
	      }
	    }
	    if(not (reported_distance <= nearest_dist)){
	      cout << i << " : " << idx << " : " << nearest << endl;
	      std::cout << "Query point: ";
	      print_point(q);
	      std::cout << "Reported neighbor: ";
	      print_point(P[idx_other]);
	      std::cout << "Reported distance: " << reported_distance << "\n";
	      std::cout << "Actual neighbor: ";
	      print_point(nearest);
	      std::cout << "Actual distance: " << nearest_dist << "\n";
	      std::cout<<"ERROR: nearest neighbor not correct"<< "\n";
	      abort();
	    }
	  }	
	}

}; //end NN_helper struct

bool do_check_correct(){
	float check = (float) rand()/RAND_MAX;
	if (check < .00001) return true;
	return false;
}

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
      return parlay::sequence<unsigned long>(k);});
  t.next("initialize answers");

  parlay::parallel_for(0, n, [&] (uint i){
      NN.ksearch(P[i], k, answers[i], eps);});
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
