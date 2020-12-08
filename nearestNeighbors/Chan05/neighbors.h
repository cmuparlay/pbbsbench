//Timothy Chan 12/05
//approximate nearest neighbors: the SSS method (static version)

#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "../octTree/oct_tree.h"
#define sq(x) (((double) (x))* ((double) (x)))
#define MAX1 (1L<<32) 



// Guy: This does not match our definition of point
typedef size_t* Point; 
size_t shift;
int key_bits = 32;
int d; 
bool report_stats = true;
bool check_correctness = true;
int algorithm_version = 0; //just because octTree/neighbors requires this parameter
double eps = 0; 

double squared_distance(Point p, Point q){
	double dist = 0; 
	for(int j = 0; j<d; j++){
		dist += sq(q[j]-p[j]);
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
void convert(parlay::sequence<vtx*> &v, size_t n, Point P[]){
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
    parlay::parallel_for(0, n, [&] (size_t i){
    	P[i] = new size_t[d];
    	for (int j = 0; j < d; j++) 
      		P[i][j] = floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta); 
    });
}

inline size_t less_msb(size_t x, size_t y) { return x < y && x < (x^y); }

//compares the interleaved bits of shifted points p, q without explicitly interleaving them
auto cmp_shuffle = [&] (Point p, Point q){
	int j, k; size_t x, y;
	for (j = k = x = 0; k<d; k++){
		if (less_msb(x, y = ((p)[k]+shift)^((q)[k]+shift))){
			j=k; 
			x=y;
		}
	}
	return (((p)[j] - (q)[j]) < 0);
};

size_t cmp_shuffle1(Point *p, Point *q){
	int j, k; size_t x, y;
	for (j = k = x = 0; k < d; k++){
		if (less_msb(x, y = ((*p)[k]+shift)^((*q)[k]+shift))){
			j = k; x = y;
		}
	}
	return (*p)[j] - (*q)[j];
}


//sort the points based on their interleave ordering
void SSS_preprocess(Point P[], size_t n){ 
	shift = 0; //TODO change back to rand
	parlay::sort_inplace(parlay::make_slice(P, P+n), cmp_shuffle);
	// std::cout << "done preprocessing" << "\n";
}

void SSS_preprocess1(Point P[], size_t n){
	shift = 0; //TODO change back to rand
	qsort((void *) P, n, sizeof(Point), 
		(int (*)(const void *, const void *)) cmp_shuffle1);
}



struct Chan_nn{

	double r, r_sq; 
	Point ans, q1, q2;  

	Chan_nn(){
		q1 = new size_t[d];
		q2 = new size_t[d];
	}

	void check_dist(Point p, Point q){
		if (are_equal(p, q)) return;
		int j; double z; 
		for (j=0, z=0; j<d; j++) z+= sq(p[j]-q[j]);
		if (z < r_sq) {
			r_sq = z; r = sqrt(z); ans = p;
			for (j=0; j<d; j++) {
				q1[j] = (q[j]>r) ? (q[j] - (size_t)ceil(r)) : 0;
				q2[j] = (q[j]+r<MAX1) ? (q[j]+(size_t)ceil(r)) : MAX1;
			}
		}
	}

	//find the distance between query point q and the box defined by p1 and p2
	double dist_sq_to_box(Point q, Point p1, Point p2){
		int i, j; size_t x, y; double z;
		for (j=x=0; j<d; j++)
			if (less_msb(x, y=(p1[j]+shift)^(p2[j]+shift))) x=y;
		frexp(x, &i);
		for (j=0, z=0; j<d; j++){
			x = ((p1[j]+shift)>>i)<<i; y = x+(1<<i);
			if (q[j]+shift < x) z += sq(q[j]+shift-x);
			else if (q[j] + shift > y) z += sq(q[j]+shift - y);
		}
		return z;
	}

	void SSS_query0(Point P[], size_t n, Point q){
		if (n==0) return;
		check_dist(P[n/2], q);
		// std::cout << "checked distance" << "\n";
		if (n == 1 || dist_sq_to_box(q, P[0], P[n-1])*sq(1+eps) > r_sq) { //getting rid of this condition did not fix the error
			return;
		}
		if (cmp_shuffle1(&q, &P[n/2]) < 0) {
			// std::cout << "cmp_shuffle result " << cmp_shuffle(q, P[n/2]) << "\n";
			SSS_query0(P, n/2, q);
			std::cout << "left "; 
			print_point(q2);
			if (cmp_shuffle1(&q2, &P[n/2]) > 0) SSS_query0(P + n/2+1, n-n/2-1, q);
		} else{
			std::cout << "right ";  
			print_point(q1);
			SSS_query0(P + n/2+1, n-n/2-1, q);
			if (cmp_shuffle1(&q1, &P[n/2]) < 0) SSS_query0(P, n/2, q);    
		}
	}

	//brute-force check correctness for approximately 100 points
	bool do_check_correct(){
		if (check_correctness){
			float check = (float) rand()/RAND_MAX;
			if (check < .00001) return true;
		}
		return false;
	}

	void check_correct(Point P[], size_t n, Point q){
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
		if(not are_equal(ans, nearest)){
			std::cout << are_equal(ans, q) << "\n";
			std::cout << "Query point: ";
			print_point(q);
			std::cout << "Reported neighbor: ";
			print_point(ans);
			std::cout << "Reported distance: " << squared_distance(q, ans) << "\n";
			std::cout << "Actual neighbor: ";
			print_point(nearest);
			std::cout << "Actual distance: " << nearest_dist << "\n";
			std::cout<<"ERROR: nearest neighbor not correct"<< "\n";
			abort();
		}		
	}

	Point SSS_query(Point P[], size_t n, Point q){
		r_sq = DBL_MAX;
		// std::cout << "reached SSS_query" << "\n";
		SSS_query0(P, n, q);
		if (not do_check_correct()) check_correct(P, n, q);		
		return ans;
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
		d = (v[0]->pt.dimension());
		Point *P;
		P = new Point[n]; 
		convert(v, n, P);
		// for (int i = 0; i<n; i++){
		// 	print_point(P[i]);
		// }
		t.next("convert to Chan's types");
		srand48(12121+n+n+d); 
		SSS_preprocess1(P, n); //parallel preprocess is buggy still
		// for (int i = 0; i<n; i++){
		// 	print_point(P[i]);
		// }
		t.next("preprocess");
		parlay::parallel_for(41, 42, [&] (size_t i){
			Chan_nn C; 
			C.SSS_query(P, n, P[i]);
		}
		);
		// for (size_t i = 0; i < n; i++){
		// 	Chan_nn C;
		// 	C.SSS_query(P, n, P[i]);
		// }
		t.next("find all");
		parlay::parallel_for(0, n, [&] (size_t i){
			delete P[i];
		});
		delete P; 
		t.next("delete data");
	};
}

