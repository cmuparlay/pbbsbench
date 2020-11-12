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
#define MAX1 (1<<29) //TODO we might need to make this bigger



// Guy: This does not match our definition of point
typedef size_t* Point; 
int d;
size_t shift;
double r, r_sq; 
Point ans, q1, q2; 
int key_bits = 64;
bool report_stats = true;
int algorithm_version = 0; //just because octTree/neighbors requires this parameter
float eps = 0; 

template<class vtx>
void convert(parlay::sequence<vtx*> &v, size_t n, int d, Point P[]){
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
    for (size_t i=0; i<n; i++){
    	P[i] = new size_t[d];
    	for (int j = 0; j < d; j++) 
      		P[i][j] = floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta); 
    }
}

inline int less_msb(size_t x, size_t y) { return x < y && x < (x^y); }

//compares the interleaved bits of shifted points p, q without explicitly interleaving them
int cmp_shuffle(Point *p, Point *q){
	int j, k; size_t x, y;
	for (j = k = x = 0; k<d; k++){
		if (less_msb(x, y = ((*p)[k]+shift)^(*q[k]+shift))){
			j=k; 
			x=y;
		}
	}
	return (*p)[j] - (*q)[k];
}


//sort the points based on their interleave ordering
void SSS_preprocess(Point P[], size_t n, int d){ 
	shift = (size_t) (drand48()*MAX1);	
	q1 = new size_t[d];
	q2 = new size_t[d];
	qsort((void *) P, n, sizeof(Point),
		(int (*)(const void *, const void *)) cmp_shuffle); //TODO can the int here remain or should it be converted to size_t?
}


void check_dist(Point p, Point q){
	int j; double z; 
	for (j=0, z=0; j<d; j++) z+= sq(p[j]-q[j]);
	if (z < r_sq) {
		r_sq = z; r = sqrt(z); ans = p;
		for (j=0; j<d; j++) {
			q1[j] = (q[j]<r) ? (q[j] - (size_t)ceil(r)): MAX1;
		}
	}
}

//find the distance between query point q and the box that p is in
float dist_sq_to_box(Point q, Point p1, Point p2){
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
	if (n == 1 || dist_sq_to_box(q, P[0], P[n-1])*sq(1+eps) > r_sq) return;
	if (cmp_shuffle(&q, &P[n/2]) > 0) {
		SSS_query0(P, n/2, q);
		if (cmp_shuffle(&q2, &P[n/2]) < 0) SSS_query0(P + n/2+1, n-n/2-1, q);
	} else{
		SSS_query0(P + n/2+1, n-n/2+1, q);
		if (cmp_shuffle(&q1, &P[n/2]) < 0) SSS_query0(P, n/2, q);
	}
}

Point SSS_query(Point P[], size_t n, Point q){
	r_sq = DBL_MAX;
	SSS_query0(P, n, q);
	return ans;
}

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
	if (k > 1) {
    	std::cout << "not implemented for k > 1" << std::endl;
    	abort();
    }
	timer t("ANN", report_stats);{
		t.next("convert to Chan's types");
		size_t n = v.size(); 
		int dims = (v[0]->pt.dimension());
		Point *P;
		P = new Point[n]; 
		convert(v, n, dims, P);
		t.next("preprocess");
		srand48(12121+n+n+dims); //TODO check this works with n as size_t
		SSS_preprocess(P, n, dims); //TODO check that this worked the way it was supposed to, maybe rounding was an issue
		t.next("find all");
		for(size_t i=0; i< n; i++){
			SSS_query(P, n, P[i]);
		}
		// parlay::parallel_for(0, n, [&] (size_t i){
		// 	SSS_query(P, n, P[i]);
		// }
		// );
	};
}

//questions about compiling using parlay:
// 1) what exactly does the makefile do here? do I just need what I have in it currently?
// 2) in testInputs, what do bnchmark and benchmark mean? why is there an identical copy in nearestNeighbors/bench?
// 3) in the octTree folder, what is the huge file just called neighbors?
// 4) it seems like some of the stuff in neighborsCheck and neighborsTime is needed here, what exactly?
