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

#include "../../octTree/oct_tree.h"

#define sq(x) (((double) (x))* ((double) (x)))

uint abs(uint x, uint y){
	if (x > y){
		return (x-y);
	} else{
		return (y-x);
	}
}

template<int d, class vtx>
struct NN_helper{
	using o_tree = oct_tree<vtx>;
	using box = typename o_tree::box;

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


	void convert(parlay::sequence<vtx*> &v, uint n, Point P[]){
	  //prelims for rounding each point to an integer: 
	  // 1) find the smallest point in each dimension
	  // 2) find the largest gap between min and max over all dimensions
	  using point = typename vtx::pointT;
	  using box = typename o_tree::box;

	  auto minmax = [&] (box x, box y) {
	    return box(x.first.minCoords(y.first),
		       x.second.maxCoords(y.second));};
	  
	  // uses a delayed sequence to avoid making a copy
	  auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
	      return box(v[i]->pt, v[i]->pt);});
	  box identity = pts[0];
	  box b = parlay::reduce(pts, parlay::make_monoid(minmax,identity));

	  double Delta = 0;
	  for (int i = 0; i < d; i++) 
	    Delta = std::max(Delta, b.second[i] - b.first[i]); 
	  point min_point = b.first;
	  parlay::sequence<vtx*> v1; 
	  v1 = parlay::sequence<vtx*>(n/2);
	  parlay::parallel_for(0, n/2, [&] (size_t i){
	  	v1[i] = v[n/2+i];
	  });
	  parlay::sequence<vtx*> v2 = z_sort(v1, b, Delta);
	  // round each point to an integer with key_bits bit
	  int bits = 31; // for some reason 32 bits does not work
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


	parlay::sequence<vtx*> z_sort(parlay::sequence<vtx*> v, box b, double Delta){ 
	    using indexed_point = typename o_tree::indexed_point; 
	    size_t n = v.size();
	    parlay::sequence<indexed_point> points;
	    points = parlay::sequence<indexed_point>(n);
	    parlay::parallel_for(0, n, [&] (size_t i){
	      size_t p1 = o_tree::interleave_bits(v[i]->pt, b.first, Delta);
	      indexed_point i1 = std::make_pair(p1, v[i]);
	      points[i] = i1; 
	    });
	    auto less = [&] (indexed_point a, indexed_point b){
	      return a.first < b.first;
	    };
	    auto x = parlay::sort(points, less);
	    parlay::sequence<vtx*> v3; 
	    v3 = parlay::sequence<vtx*>(n);
	    parlay::parallel_for(0, n, [&] (size_t i){
	      v3[i] = x[i].second; 
	    });
	    return v3; 
  	}


	void separate(uint n, Point P[], Point Q[], Point N[]){
		parlay::parallel_for(0, n, [&] (uint i){
			if((i%2)==0){
				N[i/2] = P[i/2];
			}else{
				Q[i/2] = P[i/2];    
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
