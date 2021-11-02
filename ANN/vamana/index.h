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
#include <random>
#include <set>


template<class fvec_point>
struct knn_index{
	int maxDeg;
	int beamSize;
	int k; 
	using slice_fvec = decltype(make_slice(parlay::sequence<fvec_point*>()));

	knn_index(int md, int bs, int kk) : maxDeg(md), beamSize(bs), k(kk) {}

	//give each vertex maxDeg random out neighbors
	void random_index(parlay::sequence<fvec_point*> v){
		size_t n = v.size(); 
		parlay::parallel_for(0, n, [&] (size_t i){
	    	std::random_device rd;    
			std::mt19937 rng(rd());   
			std::uniform_int_distribution<int> uni(0,n); 

			//use a set to make sure each out neighbor is unique
			std::set<int> indexset;
			while(indexset.size() < maxDeg){
				int j = uni(rng);
				indexset.insert(j);
			}
			
			for (std::set<int>::iterator it=indexset.begin(); it!=indexset.end(); ++it){
        		v[i] -> out_nbh.push_back(*it);
			} 

	    }, 1
	    );
	}


	parlay::sequence<float> centroid_helper(slice_fvec a, int d){
		if(a.size() == 1){
			parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
			for(int i=0; i<d; i++) centroid_coords[i] = (a[0]->coordinates)[i];
			return centroid_coords;
		}
		else{
			size_t n = a.size();
			// std::cout << n << std::endl; 
			parlay::sequence<float> c1;
			parlay::sequence<float> c2;
			parlay::par_do_if(n>1000,
				[&] () {c1 = centroid_helper(a.cut(0, n/2), d);},
				[&] () {c2 = centroid_helper(a.cut(n/2, n), d);}
			);
			parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
			for(int i=0; i<d; i++){
				float result = c1[i] + c2[i];
				centroid_coords[i] = result;
			}
			return centroid_coords;
		}
	}

	fvec_point* medoid_helper(fvec_point* centroid, slice_fvec a, int d){
		if(a.size() == 1){
			return a[0];
		}
		else{
			size_t n = a.size();
			fvec_point* a1;
			fvec_point* a2;
			parlay::par_do_if(n>1000,
				[&] () {a1 = medoid_helper(centroid, a.cut(0, n/2), d);},
				[&] () {a2 = medoid_helper(centroid, a.cut(n/2, n), d);}
			);
			float d1 = distance(centroid, a1);
			float d2 = distance(centroid, a2);
			if(d1<d2) return a1;
			else return a2;
		}
	}


	//computes the centroid and then assigns the approx medoid as the point in v
	//closest to the centroid
	fvec_point* find_approx_medoid(parlay::sequence<fvec_point*> v){
		size_t n = v.size();
		int d = ((v[0]->coordinates).size())/4;
		parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v), d);
		fvec_point centroidp = typename fvec_point::fvec_point();
		centroidp.coordinates = parlay::make_slice(centroid);
		fvec_point* medoid = medoid_helper(&centroidp, parlay::make_slice(v), d);
		return medoid; 
	}

	parlay::sequence<int> beam_collect(){

	}

	void build_index(parlay::sequence<fvec_point*> v){
		random_index(v);
		find_approx_medoid(v);
	}
};