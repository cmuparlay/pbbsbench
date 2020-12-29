#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
// #include "parlay/alloc.h"
#include "common/geometry.h"
#include "common/get_time.h" 
#include "../../octTree/oct_tree.h"
#include "../include/sfcnn.hpp"


int algorithm_version = 0; 
int key_bits = 64;

//this attempts to modify the sequence in-place, rounding each coordinate to a 32-bit integer
template<class vtx>
void convert(parlay::sequence<vtx*> &v, size_t n, int d, vector<uint> P[]){
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
	int bits = key_bits/d;
	int maxval = (((size_t) 1) << bits) - 1;
    parlay::parallel_for(0, n, [&] (size_t i){
    	vector<uint> vect;
    	P[i] = vect;
    	for (int j = 0; j < d; j++){
      		uint coord = (uint) floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta); 
      		vect.push_back(coord);
      	}
    });
}

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
	size_t n = v.size();
	int d = (v[0]->pt.dimension());
	vector<uint> P[n];
	convert(v, n, d, P);
	// vector<uint> *first;
	// first = P;
	sfcnn<vector<uint>, 2, uint > NN(P, n);
}