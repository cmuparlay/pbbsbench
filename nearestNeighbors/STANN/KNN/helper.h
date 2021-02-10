#define sq(x) (((double) (x))* ((double) (x)))

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
	  using point = typename vtx::pointT;
	  using box = std::pair<point,point>;

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
