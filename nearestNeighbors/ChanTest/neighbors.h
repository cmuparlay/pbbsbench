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

#define sq(x) (((float) (x))*((float) (x)))
#define MAX1 (1<<29)

//Chan's types and global variables
typedef int *Point;
int d, shift;
float eps = 0;
float r, r_sq;
Point ans, q1, q2;

//types and global variables for this particular implementation
bool report_stats = true;
int key_bits = 32;
int algorithm_version = 1;


void print_point(Point p){
	std::cout << "Point: ";
	for(int j=0; j<d; j++){
		std::cout << p[j] << ", ";
	}
	std::cout << "\n";
}


bool are_equal(Point p, Point q){
	for (int j = 0; j<d; j++){
		if(p[j] != q[j]){
			return false;
		}
	}
	return true;
}

int convert_to_int(size_t data){
	if (data > INT_MAX){
		std::cout << data << "\n";
		std::cout<<"ERROR: data point larger than INT_MAX"<< "\n";
		abort();
	}
	int convertData = static_cast<int>(data);
	return convertData;
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
    	P[i] = new int[d];
    	for (int j = 0; j < d; j++) 
      		P[i][j] = convert_to_int(floor(maxval * ((v[i] -> pt)[j] - min_point[j])/Delta)); 
    });
}



inline int less_msb(size_t x, size_t y) { return x < y && x < (x^y); }

int cmp_shuffle(Point *p, Point *q){
	int j, k, x, y;
	for (j = k = x = 0; k < d; k++){
		if (less_msb(x, y = (((*p)[k]+shift))^((*q)[k]+shift))){
			j = k; x = y;
		}
	}
	return (*p)[j]-(*q)[j];
}

void SSS_preprocess(Point P[], int n){
	shift = 0; //TODO change back to rand
	qsort((void *) P, n, sizeof(Point), 
		(int (*)(const void *, const void *)) cmp_shuffle);
}

void check_dist(Point p, Point q){
	if (are_equal(p, q)) return;
	int j; float z; 
	for (j=0, z=0; j<d; j++) z+= sq(p[j]-q[j]);
	if (z < r_sq) {
		r_sq = z; r = sqrt(z); ans = p;
		for (j=0; j<d; j++) {
			q1[j] = (q[j]>r) ? (q[j] - (int)ceil(r)) : 0;
			q2[j] = (q[j]+r<MAX1) ? (q[j]+(int)ceil(r)) : MAX1;
		}
	}
}

float dist_sq_to_box(Point q, Point p1, Point p2){
	int i, j, x, y; float z;
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

void SSS_query0(Point P[], int n, Point q){
	// std::cout << "here" << "\n";
	if (n==0) return;
	check_dist(P[n/2], q);
	// std::cout << "here" << "\n";
	if (n == 1 || dist_sq_to_box(q, P[0], P[n-1])*sq(1+eps) > r_sq) { 
		return;
	}
	if (cmp_shuffle(&q, &P[n/2]) < 0) {
		SSS_query0(P, n/2, q);
		std::cout << "left "; 
		print_point(q2);
		if (cmp_shuffle(&q2, &P[n/2]) > 0) SSS_query0(P + n/2+1, n-n/2-1, q);
	} else{
		SSS_query0(P + n/2+1, n-n/2-1, q);
		std::cout << "right "; 
		print_point(q1);
		if (cmp_shuffle(&q1, &P[n/2]) < 0) SSS_query0(P, n/2, q);
	}
}

Point SSS_query(Point P[], int n, Point q){
	r_sq = FLT_MAX;
	// std::cout << "here" << "\n";
	SSS_query0(P, n, q);
	return ans;
}

//SET OF FUNCTIONS FOR CHECKING CORRECTNESS

float squared_distance(Point p, Point q){
	float dist = 0; 
	for(int j = 0; j<d; j++){
		dist += sq(q[j]-p[j]);
	}
	return dist; 
}

void check_correct(Point P[], int n, Point q){
	Point nearest;
	float nearest_dist = FLT_MAX;
	for(int i = 0; i < n; i++){
		if(not are_equal(P[i], q)){ //make sure we don't report the query point as its own nearest neighbor
			float dist = squared_distance(q, P[i]);
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

//END CORRECTNESS CHECK 

template<int maxK, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k){
	if (k > 1) {
    	std::cout << "not implemented for k > 1" << std::endl;
    	abort();
    }
	timer t("ANN", report_stats);{
		size_t n = v.size(); 
		int m = convert_to_int(n);
		d = (v[0]->pt.dimension());
		Point *P;
		P = new Point[m]; 
		convert(v, n, P);
		// for (int i=0; i< m; i++){
		// 	print_point(P[i]);
		// }
		t.next("convert to Chan's types");
		srand48(12121+m+m+d); 
		SSS_preprocess(P, m); 
		t.next("preprocess");
		// for (int i=0; i< m; i++){
		// 	print_point(P[i]);
		// }
		q1 = new int[d];
		q2 = new int[d];
		for (int i =0; i<1; i++){
			SSS_query(P, m, P[i]);
			check_correct(P, m, P[i]);	
		}
		t.next("find all");
		for (int i = 0; i<m; i++) delete P[i];
		delete P; 
		t.next("delete data");
	};
}