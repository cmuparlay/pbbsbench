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
#include "parlay/random.h"
#include "../utils/indexTools.h"
#include "common/geometry.h"
#include <random>
#include <set>
#include <math.h>

extern bool report_stats;


template<typename T>
struct knn_index{
	int maxDeg; // maximum out-degree of a vertex in the graph
	int beamSize; // beam size for beam search during querying and construction
	double r2_alpha; // alpha parameter for round 2 of robustPrune (or if there's no round 2, for round 1)
	unsigned d; // the dimension of the points
	bool mips; // whether the search is maximum inner product search as opposed to minimal euclidean distance
	std::set<int> delete_set; // for use in lazy deletion of points
	using tvec_point = Tvec_point<T>;
	using fvec_point = Tvec_point<float>;
	tvec_point* medoid;
	using pid = std::pair<int, float>;
	using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
	using index_pair = std::pair<int, int>;
	using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));

	knn_index(int md, int bs, double a, unsigned dim, bool m=false) : maxDeg(md), beamSize(bs), r2_alpha(a), d(dim), mips(m) {}

	/* 
	function for computing the distance between two points in the dataset

	@param p: pointer to the first point
	@param q: pointer to the second point
	@param d: dimension of the points
	*/
	float Distance(T* p, T* q, unsigned d){
		if(mips) return mips_distance(p, q, d);
		else return distance(p, q, d);
	}

	/* 
	function for computing the distance between a query and an element of the dataset 
	
	@param p: pointer to the query
	@param q: pointer to the point to be compared to the query
	@param d: dimension of the points
	*/
	float CDistance(float* p, T* q, unsigned d){
		if(mips) return mips_distance(p, q, d);
		else return distance(p, q, d);
	}

	parlay::sequence<float> centroid_helper(slice_tvec a){
		if(a.size() == 1){
			parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
			for(int i=0; i<d; i++) centroid_coords[i] = static_cast<float>((a[0]->coordinates)[i]);
			return centroid_coords;
		}
		else{
			size_t n = a.size();
			parlay::sequence<float> c1;
			parlay::sequence<float> c2;
			parlay::par_do_if(n>1000,
				[&] () {c1 = centroid_helper(a.cut(0, n/2));},
				[&] () {c2 = centroid_helper(a.cut(n/2, n));}
			);
			parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
			for(int i=0; i<d; i++){
				float result = (c1[i] + c2[i])/2;
				centroid_coords[i] = result;
			}
			return centroid_coords;
		}
	}

	tvec_point* medoid_helper(fvec_point* centroid, slice_tvec a){
		if(a.size() == 1){
			return a[0];
		}
		else{
			size_t n = a.size();
			tvec_point* a1;
			tvec_point* a2;
			parlay::par_do_if(n>1000,
				[&] () {a1 = medoid_helper(centroid, a.cut(0, n/2));},
				[&] () {a2 = medoid_helper(centroid, a.cut(n/2, n));}
			);
			float d1 = CDistance(centroid->coordinates.begin(), a1->coordinates.begin(), d);
			float d2 = CDistance(centroid->coordinates.begin(), a2->coordinates.begin(), d);
			if(d1<d2) return a1;
			else return a2;
		}
	}

	//computes the centroid and then assigns the approx medoid as the point in v
	//closest to the centroid
	void find_approx_medoid(parlay::sequence<Tvec_point<T>*> &v){
		size_t n = v.size();
		parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
		fvec_point centroidp = Tvec_point<float>();
		centroidp.coordinates = parlay::make_slice(centroid);
		medoid = medoid_helper(&centroidp, parlay::make_slice(v));
		std::cout << "Medoid ID: " << medoid->id << std::endl;
	}

	int get_medoid(){return medoid->id;}

	void print_set(std::set<int> myset){
		std::cout << "[";
		for (std::set<int>::iterator it=myset.begin(); it!=myset.end(); ++it){
			std::cout << *it << ", ";
		}
		std::cout << "]" << std::endl;
	}

	
	/* 
	robustPrune routine as found in DiskANN paper, with the exception that the new candidate set
	is added to the field new_nbhs instead of directly replacing the out_nbh of p to allow concurrent pruning without locks.

	@param p: the point whose out_nbh is being pruned/updated
	@param candidates: the candidate set, consisting of pid pairs with indices of vectors and their distance to p (per the paper, all points visited in beam search from medoid to p over the existing graph). Not expected to contain the neighbors of p iff add is false. If add is false and the neighbors are not contained in the candidate set, all current neighbors of p will be removed (or at least the ones not incidentally in the candidate set).
	@param v: the dataset, a sequence containing pointers to the tvec_point objects which represent the vectors/graph vertices
	@param alpha: the pruning parameter
	@param add: whether or not to add the out neighbors of p to the candidate set (see candidates)
	*/
	void robustPrune(tvec_point* p, parlay::sequence<pid> candidates, parlay::sequence<tvec_point*> &v, double alpha, bool add = true) {
    // add out neighbors of p to the candidate set.
    if(add){
    	for (size_t i=0; i < size_of(p->out_nbh); i++) {
				candidates.push_back(std::make_pair(p->out_nbh[i],
					Distance(v[p->out_nbh[i]]->coordinates.begin(), p->coordinates.begin(), d)));
			}
    }
		
		// Sort the candidate set in reverse order according to distance from p.
		auto less = [&] (pid a, pid b) {return a.second < b.second;};
		parlay::sort_inplace(candidates, less);

		parlay::sequence<int> new_nbhs = parlay::sequence<int>(); // new neighbors of p
    // new_nbhs.reserve(maxDeg);

    size_t candidate_idx = 0; // index of the candidate to be considered in the sorted sequence
		while (new_nbhs.size() < maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      int p_star = candidates[candidate_idx].first; // the id of the candidate to be considered
      candidate_idx++;
      if (p_star == p->id || p_star == -1) { // if the candidate is p itself or has been pruned
        continue;
      }

      new_nbhs.push_back(p_star); // add current best candidate to new_nbhs
	  // Prune the candidate set.
      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        int p_prime = candidates[i].first;
        if (p_prime != -1) {
          float dist_starprime = Distance(v[p_star]->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = -1; // pruned candidates have -1 as their id in the candidate set
          }
        }
      }
		}
		add_new_nbh(new_nbhs, p);
	}


	/* 
	wrapper to allow calling robustPrune on a sequence of candidates that do not come with precomputed distances 

	@param p: the point whose out_nbh is being pruned/updated
	@param candidates: the candidate set, consisting of indices of vectors (per the paper, all points visited in beam search from medoid to p over the existing graph). Not expected to contain the neighbors of p iff add is false. If add is false and the neighbors are not contained in the candidate set, all current neighbors of p will be removed (or at least the ones not incidentally in the candidate set).
	@param v: the sequence containing all points already in the graph
	@param alpha: the pruning parameter
	*/
	void robustPrune(Tvec_point<T>* p, parlay::sequence<int> candidates, parlay::sequence<Tvec_point<T>*> &v, double alpha, bool add = true){
    parlay::sequence<pid> cc;
    cc.reserve(candidates.size() + size_of(p->out_nbh));
    for (size_t i=0; i<candidates.size(); ++i) {
      cc.push_back(std::make_pair(candidates[i], Distance(v[candidates[i]]->coordinates.begin(), p->coordinates.begin(), d)));
    }
    return robustPrune(p, std::move(cc), v, alpha, add);
	}

	/*
	constructs the graph with dataset v, with inserts describing the elements of the dataset to be inserted into the graph.
	*/
	void build_index(parlay::sequence<Tvec_point<T>*> &v, parlay::sequence<int> inserts, bool two_pass=false){
		std::cout << "Mips: " << mips << std::endl;
		clear(v);
		find_approx_medoid(v);
		if(two_pass){
		  std::cout << "Starting first pass" << std::endl; 
		  batch_insert(inserts, v, true, 1.0, 2, .02, two_pass);
		  std::cout << "Finished first pass" << std::endl;
		} 
		batch_insert(inserts, v, true, r2_alpha, 2, .02,  two_pass);
	}

	/* 
	Adds a sequence of indices representing points to the delete_set, which will exclude them from query results until consolidate_deletes() is called and they're actually removed from the graph
	*/
	void lazy_delete(parlay::sequence<int> deletes, parlay::sequence<Tvec_point<T>*> &v){
		for(int p : deletes){
			if(p < 0 || p > (int) v.size() ){
				std::cout << "ERROR: invalid point " << p << " given to lazy_delete" << std::endl; 
				abort();
			}
			if(p != medoid->id) delete_set.insert(p);
			else std::cout << "Deleting medoid not permitted; continuing" << std::endl; 
		} 
	}

	/*
	adds an index representing a point to the delete_set, which will exclude it from query results until consolidate_deletes() is called and it's actually removed from the graph
	*/
	void lazy_delete(int p, parlay::sequence<Tvec_point<T>*> &v){
		if(p < 0 || p > (int) v.size()){
			std::cout << "ERROR: invalid point " << p << " given to lazy_delete" << std::endl; 
			abort();
		}
		if(p == medoid->id){
			std::cout << "Deleting medoid not permitted; continuing" << std::endl; 
			return;
		} 
		delete_set.insert(p);
	}

	/* 
	The Delete algorithm described in Algorithm 4 of the FreshDiskANN paper. Updates the graph to remove all points in the delete_set (while preserving the alpha-RNG property), and then clears the delete_set.
	*/
	void consolidate_deletes(parlay::sequence<Tvec_point<T>*> &v){
		//clear deleted neighbors out of delete set for preprocessing
		parlay::parallel_for(0, v.size(), [&] (size_t i){ // for each point i in the graph
			if(delete_set.find(i) != delete_set.end()){ // if i is in the delete set
				parlay::sequence<int> new_edges; 
				for(int j=0; j<size_of(v[i]->out_nbh); j++){ // for each neighbor of i
					if(delete_set.find(v[i]->out_nbh[j]) == delete_set.end()) new_edges.push_back(v[i]->out_nbh[j]); // if said neighbor is not in the delete set, add it to new_edges (the new out_nbh); don't delete it
				}
				if(new_edges.size() < size_of(v[i]->out_nbh)) add_out_nbh(new_edges, v[i]); // if the neighborhood has changed (neighbors were deleted), update it
			}
		});

		parlay::parallel_for(0, v.size(), [&] (size_t i){ // for each point i in the graph (again)
			if(delete_set.find(i) == delete_set.end() && size_of(v[i]->out_nbh) != 0){ // if i is not in the delete set and has a non-empty out_nbh
				std::set<int> new_edges;
				bool modify = false;
				for(int j=0; j<size_of(v[i]->out_nbh); j++){ // for each neighbor j of i
					int index = v[i]->out_nbh[j]; // (moved this from the else -Ben)
					if(delete_set.find(index) == delete_set.end()){ 
						new_edges.insert(index); // if j is not in the delete set add j to new_edges
					} else{ 
						modify = true; // now that we're deleting a neighbor the neighborhood will need to be replaced by new_edges
						for(int k=0; k<size_of(v[index]->out_nbh); k++) new_edges.insert(v[index]->out_nbh[k]); // add all of j's neighbors to new_edges
					}
				}
				//TODO only prune if overflow happens
				//TODO modify in separate step with new memory initialized in one block
				if(modify){ 
					parlay::sequence<int> candidates; 
					for(int j : new_edges) candidates.push_back(j); // add all new_edges to candidates (this is essentially just converting the new_edges set to a sequence, but it's not obvious why the set would be faster than a preallocated sequence (size_of(v[i]->out_nbh)) to begin with if this is the only time the values are accessed)
					parlay::sequence<int> new_neighbors(maxDeg, -1); // constructs a sequence of size maxDeg with all elements -1
					v[i]->new_nbh = parlay::make_slice(new_neighbors.begin(), new_neighbors.end());
					robustPrune(v[i], candidates, v, r2_alpha, false); // prune candidates to get i's new neighbors
					synchronize(v[i]);
				}
				
			} 
		});
		// for every point i in the graph, if it's in the delete set set all its out neighbors to -1
		// is it actually faster to iterate in parallel over the whole graph than to iterate over the delete set?
		parlay::parallel_for(0, v.size(), [&] (size_t i){
			if(delete_set.find(i) != delete_set.end()){
				clear(v[i]);
			} 
		});
 
		delete_set.clear(); // clear the delete set

	}

	void insert_and_count(parlay::sequence<int> &inserts, parlay::sequence<Tvec_point<T>*> &v, 
		parlay::sequence<std::pair<int, int>> &edges_to_record){

		double alpha = 1.2;

		size_t floor = 0; // I don't think either of these are actually reassigned within this function, which I think makes at least floor redundant
		size_t ceiling = inserts.size(); 
		
		// does this actually shuffle the sequence? This seems like a good use for parlay::random_permutation
		auto shuffled_inserts = parlay::tabulate(inserts.size(), [&] (size_t i){return inserts[i];});

		parlay::sequence<int> new_out = parlay::sequence<int>(maxDeg*(ceiling-floor), -1); // this is a sequence of size maxDeg*(ceiling-floor) with all elements -1
		//search for each node starting from the medoid, then call
		//robustPrune with the visited list as its candidate set
		parlay::parallel_for(floor, ceiling, [&] (size_t i){
			size_t index = shuffled_inserts[i];
			v[index]->new_nbh = parlay::make_slice(new_out.begin()+maxDeg*(i-floor), new_out.begin()+maxDeg*(i+1-floor)); // initializing new_nbh to a slice of new_out, which are for all intents and purposes random vertices
			parlay::sequence<pid> visited = (beam_search(v[index], v, medoid, beamSize, d, mips)).first.second; // this signature of beam_search returns a pair with a pair of sequences as first and an int as second, with the pair of sequences being the candidate list and visited list, so this ultimately returns the visited list
			if(report_stats) v[index]->visited = visited.size();
			robustPrune(v[index], visited, v, alpha);
		});

		//make each edge bidirectional by first adding each new edge
		//(i,j) to a sequence, then semisorting the sequence by key values
		auto to_flatten = parlay::tabulate(ceiling-floor, [&] (size_t i){
			int index = shuffled_inserts[i+floor];
			auto edges = parlay::tabulate(size_of(v[index]->new_nbh), [&] (size_t j){
				return std::make_pair((v[index]->new_nbh)[j], index);
			}); // is there a performance penalty for nesting these parlay functions when the length of the outer tabulate is >> p? There are earlier places where it seems these nested parallel calls are avoided.
			return edges;
		});

		auto e_1 = parlay::tabulate(ceiling-floor, [&] (size_t i) {
			int index = shuffled_inserts[i+floor];
			auto edges = parlay::tabulate(size_of(v[index]->new_nbh), [&] (size_t j){
				return std::make_pair((v[index]->new_nbh)[j], 1);
			});
			return edges;
		});
		edges_to_record = parlay::flatten(e_1);
		parlay::parallel_for(floor, ceiling, [&] (size_t i) {synchronize(v[shuffled_inserts[i]]);} );
		auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
		//finally, add the bidirectional edges; if they do not make
		//the vertex exceed the degree bound, just add them to out_nbhs;
		//otherwise, use robustPrune on the vertex with user-specified alpha
		parlay::parallel_for(0, grouped_by.size(), [&] (size_t j){
			auto [index, candidates] = grouped_by[j];
			int newsize = candidates.size() + size_of(v[index]->out_nbh);
			if(newsize <= maxDeg){
				for(const int& k : candidates) add_nbh(k, v[index]);
			} else{
				parlay::sequence<int> new_out_2(maxDeg, -1);
				v[index]->new_nbh=parlay::make_slice(new_out_2.begin(), new_out_2.begin()+maxDeg);
				robustPrune(v[index], candidates, v, r2_alpha);  
				synchronize(v[index]);
			}
		});	
	}
	/* 
	Inserts a set of points into the graph in batches of exponentially increasing size, allowing for an efficient lock-less parallel construction as a batch of points only considers the points that are already in the graph (1/base of all the points in the graph up to that point) while constructing its own out edges.

	@param inserts: a sequence of indices of points to insert into the graph
	@param v: a sequence of pointers to the points of the dataset/ vertices of the graph
	@param random_order: if true, the points are inserted in a random order
	@param alpha: the alpha parameter for robustPrune
	@param base: the base of the exponential batch size increase
	@param max_fraction: the maximum fraction of the points in the full dataset that can be inserted in a single batch
	@param two_pass: whether this is the second pass of a two-pass construction
	*/
	void batch_insert(parlay::sequence<int> &inserts, parlay::sequence<Tvec_point<T>*> &v, bool random_order = false, double alpha = 1.2, double base = 2,
		double max_fraction = .02, bool two_pass = false){
		std::cout << "alpha " << alpha << std::endl; 
		if(two_pass == false){
			for(int p : inserts){
				if(p < 0 || p > (int) v.size() || (v[p]->out_nbh[0] != -1 && v[p]->id != medoid->id)){
					std::cout << "ERROR: invalid or already inserted point " << p << " given to batch_insert" << std::endl;
					abort();
				}
			}
		}
		size_t n = v.size();
		size_t m = inserts.size();
		size_t inc = 0;
		size_t count = 0;
		size_t max_batch_size = std::min(static_cast<size_t>(max_fraction*static_cast<float>(n)), 1000000ul);
		parlay::sequence<int> rperm;
		if(random_order) rperm = parlay::random_permutation<int>(static_cast<int>(m));
		else rperm = parlay::tabulate(m, [&] (int i) {return i;});
		auto shuffled_inserts = parlay::tabulate(m, [&] (size_t i) {return inserts[rperm[i]];});
		while(count < m){
			size_t floor;
			size_t ceiling;
			if(pow(base,inc) <= max_batch_size){ // batch size increases exponentially 
				floor = static_cast<size_t>(pow(base, inc)) - 1;
				ceiling = std::min(static_cast<size_t>(pow(base, inc+1)), m) - 1;
				count = std::min(static_cast<size_t>(pow(base, inc+1)), m) - 1;
			} else{
				floor = count;
				ceiling = std::min(count + static_cast<size_t>(max_batch_size), m) - 1;
				count += static_cast<size_t>(max_batch_size);
			}
			parlay::sequence<int> new_out = parlay::sequence<int>(maxDeg*(ceiling-floor), -1);
			//search for each node starting from the medoid, then call
			//robustPrune with the visited list as its candidate set
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				size_t index = shuffled_inserts[i];
				v[index]->new_nbh = parlay::make_slice(new_out.begin()+maxDeg*(i-floor), new_out.begin()+maxDeg*(i+1-floor));
				parlay::sequence<pid> visited = (beam_search(v[index], v, medoid, beamSize, d, mips)).first.second;
				if(report_stats) v[index]->visited = visited.size();
				robustPrune(v[index], visited, v, alpha);
			});
			//make each edge bidirectional by first adding each new edge
			//(i,j) to a sequence, then semisorting the sequence by key values
			auto to_flatten = parlay::tabulate(ceiling-floor, [&] (size_t i){
				int index = shuffled_inserts[i+floor];
				auto edges = parlay::tabulate(size_of(v[index]->new_nbh), [&] (size_t j){
					return std::make_pair((v[index]->new_nbh)[j], index);
				});
				return edges;
			});
			// below, the new_nbh field replaces the old neighbors of the vertex, applying the updates determined previously
			parlay::parallel_for(floor, ceiling, [&] (size_t i) {synchronize(v[shuffled_inserts[i]]);} );
			auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
			//finally, add the bidirectional edges; if they do not make
			//the vertex exceed the degree bound, just add them to out_nbhs;
			//otherwise, use robustPrune on the vertex with user-specified alpha
			parlay::parallel_for(0, grouped_by.size(), [&] (size_t j){
				auto [index, candidates] = grouped_by[j];
				int newsize = candidates.size() + size_of(v[index]->out_nbh);
				if(newsize <= maxDeg){
					for(const int& k : candidates) add_nbh(k, v[index]);
				} else{
					parlay::sequence<int> new_out_2(maxDeg, -1);
					v[index]->new_nbh=parlay::make_slice(new_out_2.begin(), new_out_2.begin()+maxDeg);
					robustPrune(v[index], candidates, v, r2_alpha);  
					synchronize(v[index]);
				}
			});
			inc += 1;
		}
	}

	/* 
	Inserts a single point into the graph
	*/
	void batch_insert(int p, parlay::sequence<Tvec_point<T>*> &v){
		parlay::sequence<int> inserts;
		inserts.push_back(p);
		batch_insert(inserts, v, true);
	}

	/* 
	Confirms that a point has been deleted from the graph by checking that it is not itself in the graph or a neighbor of any other point
	*/
	void check_index(parlay::sequence<Tvec_point<T>*> &v){
		parlay::parallel_for(0, v.size(), [&] (size_t i){
      if(v[i]->id > 1000000 && v[i]->id != medoid->id){
      	if(size_of(v[i]->out_nbh) != 0) {
      		std::cout << "ERROR : deleted point " << i << " still in graph" << std::endl; 
      		abort();
      	}
      }else{
      	for(int j=0; j<size_of(v[i]->out_nbh); j++){
      		int nbh = v[i]->out_nbh[j];
      		if(nbh > 1000000 && nbh != medoid->id){
      			std::cout << "ERROR : point " << i << " contains deleted neighbor " << nbh << std::endl; 
      			abort();
      		}
      	}
      }
    });
	}


  void searchNeighbors(parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<Tvec_point<T>*> &v, int beamSizeQ, int k, float cut){
    searchAll(q, v, beamSizeQ, k, d, medoid, mips, cut);
  }

  void rangeSearch(parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<Tvec_point<T>*> &v, 
  	int beamSizeQ, double r, int k, float cut=1.14, double slack = 3.0){
		rangeSearchAll(q, v, beamSizeQ, d, medoid, r, k, cut, slack);
	}
};
