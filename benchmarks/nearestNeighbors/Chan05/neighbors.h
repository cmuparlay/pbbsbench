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

//This code is based on Timothy Chan's 2006 paper,  
//"A minimalist's implementation of an approximate nearest neighbor algorithm in fixed dimensions" 

#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/alloc.h"
#include "common/geometry.h"
#include "../octTree/oct_tree.h"
#include "counter.h"
#define sq(x) (((double) (x))* ((double) (x)))

typedef size_t* Point; 
typedef double* DPoint;
size_t shift, MAX1;
int key_bits = 60;
int d; 
bool report_stats = true;
bool check_correctness = false; //optional correctness check; very slow
int algorithm_version = 0; //just to comply with the makefile requirements
double eps = 0; 
bool midpoint_search = true;
int node_cutoff = 32;
bool leaf_search = false;

//because std::abs doesn't work on size_t type
size_t abs(size_t x, size_t y){
	if (x > y){
		return (x-y);
	} else{
		return (y-x);
	}
}

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
std::pair<DPoint, double> convert(parlay::sequence<vtx*> &v, size_t n, Point P[]){
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
	// add each integer point to Chan's array
	int bits = key_bits/d;
	size_t maxval = (((size_t) 1) << bits) - 1;
	MAX1 = maxval; 
    parlay::parallel_for(0, n, [&] (size_t i){
    	P[i] = (Point) parlay::p_malloc(65);
	    	for (int j = 0; j < d; j++) 
	      		P[i][j] = (size_t) floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta);  	
    });
	DPoint min_Point;
	min_Point = (DPoint) parlay::p_malloc(65*d);
	for(int i=0; i<d; i++) {
		min_Point[i] = min_point[i];
		// std::cout << min_Point[i] << std::endl;
	}
	return std::make_pair(min_Point, Delta);
}

inline bool less_msb(size_t x, size_t y) { return x < y && x < (x^y); }

//compares the interleaved bits of shifted points p, q without explicitly interleaving them
auto cmp_shuffle = [] (Point p, Point q){
	int j, k; size_t x, y;
	for (j = k = x = 0; k<d; k++){
		if (less_msb(x, y = ((p)[k]+shift)^((q)[k]+shift))){
			j=k; 
			x=y;
		}
	}
	return (((p)[j] < (q)[j]));     
};          

//sort the points based on their interleave ordering
void SSS_preprocess(Point P[], size_t n){ 
	shift = 0; //can be changed to 0 (in theory this voids the time guarantees but no difference in practice)
	parlay::sort_inplace(parlay::make_slice(P, P+n), cmp_shuffle);}



template<int d, class vtx>
struct Chan_nn{

	double r, r_sq; 
	Point ans;

	DPoint min_point;
	double Delta;

	size_t q1_[d], q2_[d];
	Point q1 = q1_;
	Point q2 = q2_; 
	
	atomic_sum_counter<size_t> recursive_calls;

	parlay::sequence<vtx*> v;

	Chan_nn(DPoint min, double D){
		min_point = min;
		Delta = D;
		r_sq = DBL_MAX;
		recursive_calls.reset();
	}

	void check_dist(Point p, Point q){
		if (are_equal(p, q)) return;
		int j; double z; 
		for (j=0, z=0; j<d; j++) z+= sq(abs(p[j], q[j])); 
		if (z < r_sq) {
			r_sq = z; r = sqrt(z); ans = p;
			for (j=0; j<d; j++) {
				if(q[j]>r){
					if(q[j] < (size_t)ceil(r)){ 
						q1[j] = 0;
					} else q1[j] = (q[j]-(size_t)ceil(r));
				} else q1[j] = 0;
				q2[j] = (q[j]+r<MAX1) ? (q[j]+(size_t)ceil(r)) : MAX1; //MAX is dependent on key_bits now
			}
		}
	}



  using point = typename vtx::pointT;
  using box = std::pair<point, point>;
  using indexed_point = std::pair<size_t, vtx*>;


  template <typename Seq>
  static box get_box(Seq &V) { 
    size_t n = V.size();
    auto minmax = [&] (box x, box y) {
      return box(x.first.minCoords(y.first),
     x.second.maxCoords(y.second));};
    // uses a delayed sequence to avoid making a copy
    auto pts = parlay::delayed_seq<box>(n, [&] (size_t i) {
      return box(V[i]->pt, V[i]->pt);});
    box identity = pts[0];  
    box final = parlay::reduce(pts, parlay::make_monoid(minmax, identity));
    return (final);
  }


	

	//find the distance between query point q and the box defined by p1 and p2
	double dist_sq_to_box(Point q, Point p1, Point p2){
		int i, j; size_t x, y; double z;
		for (j=x=0; j<d; j++)
			if (less_msb(x, y=(p1[j]+shift)^(p2[j]+shift))) x=y;
		frexp(x, &i);
		for (j=0, z=0; j<d; j++){
			x = ((p1[j]+shift)>>i)<<i; y = x+(1<<i);
			if (q[j]+shift < x) z += sq(abs(q[j] + shift, x)); 
			else if (q[j] + shift > y) z += sq(q[j]+shift - y); //this is safe because we know the result will be positive
		}
		return z;
	}

	void SSS_query0(Point P[], size_t l, size_t h, Point q){
		size_t m = h-l;
		size_t mid = l+m/2;
		if (m==0) return;
		recursive_calls.update_value(1);
		check_dist(P[mid], q);
		if (m == 1 || dist_sq_to_box(q, P[l], P[h-1])*sq(1+eps) > r_sq) { 
			return;
		} 
		if(leaf_search && m < node_cutoff) brute_force(P, l, h, q);
		if(cmp_shuffle(q, P[mid])){
			SSS_query0(P, l, mid, q);
			if(not cmp_shuffle (q2, P[mid])) SSS_query0(P, mid+1, h, q);
		} else{
			SSS_query0(P, mid+1, h, q);
			if(cmp_shuffle(q1, P[mid])) SSS_query0(P, l, mid, q);
		}
	}

	void brute_force(Point P[], size_t l, size_t h, Point q){
		for(size_t i=l; i<h; i++) check_dist(P[i], q);
	}

	void SSS_query_midchoice(Point P[], size_t l, size_t h, Point q, 
		parlay::slice<indexed_point*, indexed_point*> tags, int bit){
		if(h-l == 0) return;
		if(bit == 0) return brute_force(P, l, h, q);
		if(h-l < node_cutoff && leaf_search) return brute_force(P, l, h, q);
		else{
			// std::cout << "here1" << std::endl;
			// std::cout << bit << std::endl;
			size_t val = ((size_t) 1) << (bit - 1);
			size_t mask = (bit == 64) ? ~((size_t) 0) : ~(~((size_t) 0) << bit);
			auto less = [&] (indexed_point x) { return (x.first & mask) < val; };
			// std::cout << tags.size() << std::endl;
			// std::cout << "binary searching" << std::endl;
			// and then we binary search for the cut point
			size_t pos = parlay::internal::binary_search(tags, less);
			// std::cout << pos << std::endl;
			// std::cout << "here2" << std::endl;
			if(pos == 0 || pos == h-l) return SSS_query_midchoice(P, l, h, q, tags, bit-1);
			size_t split = l+pos-1;
			// std::cout << "here3" << std::endl;
			check_dist(P[split], q);
			recursive_calls.update_value(1);
			// std::cout << "here4" << std::endl;
			if(h-l == 1 || dist_sq_to_box(q, P[l], P[h-1])*sq(1+eps) > r_sq) return;
			if(cmp_shuffle(q, P[split])){
				// std::cout << "here5" << std::endl;
				auto l_tags = tags.cut(0, pos);
				// std::cout << "size of ltags " << l_tags.size() << std::endl;
				SSS_query_midchoice(P, l, split+1, q, l_tags, bit-1);
				if(not cmp_shuffle (q2, P[split])){ 
					auto r_tags = tags.cut(pos, h-l);
					// std::cout << "size of rtags " << r_tags.size() << std::endl;
					SSS_query_midchoice(P, split+1, h, q, r_tags, bit-1);
				}
			} else{
				// std::cout << "here6" << std::endl;
				auto r_tags = tags.cut(pos, h-l);
				// std::cout << pos << std::endl;
				// std::cout << "size of rtags " << r_tags.size() << std::endl;
				SSS_query_midchoice(P, split+1, h, q, r_tags, bit-1);
				if(cmp_shuffle(q1, P[split])){
					auto l_tags = tags.cut(0, pos);
					// std::cout << "size of ltags " << l_tags.size() << std::endl;
					SSS_query_midchoice(P, l, split+1, q, l_tags, bit-1);
				}
			}
		}
		

		
	}


	//brute-force check correctness for approximately 100 points out of every 10 million
	bool do_check_correct(){
		float check = (float) rand()/RAND_MAX;
		if (check < .1) return true;
		return false;
	}

	void check_correct(Point P[], size_t n, Point q){
		if (do_check_correct()){
			Point nearest;
			double nearest_dist = DBL_MAX;
			for(size_t i = 0; i < n; i++){
				if(not are_equal(P[i], q)){ //make sure we don't report the query point as its own nearest neighbor
					double dist = squared_distance(q, P[i]);
					if(dist < nearest_dist){
						nearest = P[i];
						nearest_dist = dist; 
					}
				}
			} 
			if(not (r_sq <= (1+eps)*nearest_dist)){
				std::cout << "Query point: ";
				print_point(q);
				std::cout << "Reported neighbor: ";
				print_point(ans);
				std::cout << "Reported distance: " << squared_distance(q, ans) << "\n";
				std::cout << "Actual neighbor: ";
				print_point(nearest);
				std::cout << "Actual distance: " << (1+eps)*nearest_dist << "\n";
				std::cout<<"ERROR: nearest neighbor not correct"<< "\n";
				abort();
			}
		}		
	}

	Point SSS_query(Point P[], size_t n, Point q, parlay::sequence<indexed_point> &tags){
		// parlay::slice<indexed_point*, indexed_point*> tagged = parlay::make_slice(tags);
		if(midpoint_search) SSS_query0(P, 0, n, q);
		else{
			// std::cout << "making slice" << std::endl;
			parlay::slice<indexed_point*, indexed_point*> tagged = parlay::make_slice(tags);
			// std::cout << "entering recursion" << std::endl;
			SSS_query_midchoice(P, 0, n, q, tagged, key_bits);
		}
		return ans;
	}

	size_t num_calls(){
		return recursive_calls.get_value();
	}

}; //end Chan_nn struct



template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
	if (k > 1) {
    	std::cout << "not implemented for k > 1" << std::endl;
    	abort();
    }
	timer t("ANN", report_stats);{
		size_t n = v.size(); 
		size_t num_trials = n;
		d = (v[0]->pt.dimension());
		parlay::sequence<size_t> recursive_calls(n, 0);
		Point *P;
		P = (Point*) parlay::p_malloc((n+1)*64); 
		auto [min_point, Delta] = convert(v, n, P);
		// for(int i=0; i<d; i++) std::cout << min_point[i] << std::endl;
		using o_tree = oct_tree<vtx>;
		// auto tagged_points = o_tree::tag_points_external(v);
		using indexed_point = std::pair<size_t, vtx*>;
		parlay::sequence<indexed_point> tagged_points;
		t.next("convert to Chan's types");
		srand48(12121+n+n+d); 
		SSS_preprocess(P, n); 
		t.next("preprocess");
		parlay::parallel_for(0, num_trials, [&] (size_t i){
			if(d==2){
				// std::cout << "NEW QUERY" << std::endl;
				// Chan_nn<2, vtx> C;
				using CNN = Chan_nn<2, vtx>;
				CNN C(min_point, Delta);
				C.SSS_query(P, n, P[i], tagged_points);
				recursive_calls[i] = C.num_calls();
			} else{
				// Chan_nn<3, vtx> C; 
				using CNN = Chan_nn<3, vtx>;
				CNN C(min_point, Delta);
				C.SSS_query(P, n, P[i], tagged_points);
				recursive_calls[i] = C.num_calls();
			}
		}
		);
		t.next("find all");
		if (check_correctness){
			parlay::parallel_for(0, num_trials, [&] (size_t i){
				if(d==2){
					// Chan_nn<2, vtx> C;
					using CNN = Chan_nn<2, vtx>;
					CNN C(min_point, Delta);
					C.SSS_query(P, n, P[i], tagged_points);
					C.check_correct(P, n, P[i]);
					recursive_calls[i] = C.num_calls();
				} else{
					// Chan_nn<3, vtx> C;
					using CNN = Chan_nn<3, vtx>;
					CNN C(min_point, Delta); 
					C.SSS_query(P, n, P[i], tagged_points);
					C.check_correct(P, n, P[i]);
					recursive_calls[i] = C.num_calls();
				}
			}
			);
			t.next("check correctness");
		}
		size_t total = parlay::reduce(recursive_calls);
		std::cout << "Total recursive calls " << total << std::endl;
		std::cout << "Avg recursive calls " << ( (float) total/ (float) num_trials) << std::endl;

		parlay::parallel_for(0, n, [&] (size_t i){
			parlay::p_free(P[i]);
		});
		parlay::p_free(P); 
		t.next("delete data");
	};
}

