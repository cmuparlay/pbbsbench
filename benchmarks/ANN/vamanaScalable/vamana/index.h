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
#include "../utils/seq_allocator.h"
#include "../utils/indexTools.h"
#include "common/geometry.h"
#include <random>
#include <set>
#include <math.h>

extern bool report_stats;


template<typename T>
struct knn_index{
	int maxDeg;
	int beamSize;
	// int k;
	double r2_alpha; //alpha parameter for round 2 of robustPrune
	unsigned d;
	using tvec_point = Tvec_point<T>;
	using fvec_point = Tvec_point<float>;
	tvec_point* medoid;
	using pid = std::pair<int, float>;
	using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
	using index_pair = std::pair<int, int>;
	using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));
	using fine_sequence = parlay::sequence<int>;

	knn_index(int md, int bs, double a, unsigned dim) : maxDeg(md), beamSize(bs), r2_alpha(a), d(dim) {}

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
			float d1 = distance(centroid->coordinates.begin(), a1->coordinates.begin(), d);
			float d2 = distance(centroid->coordinates.begin(), a2->coordinates.begin(), d);
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

	void print_set(std::set<int> myset){
		std::cout << "[";
		for (std::set<int>::iterator it=myset.begin(); it!=myset.end(); ++it){
			std::cout << *it << ", ";
		}
		std::cout << "]" << std::endl;
	}

	//robustPrune routine as found in DiskANN paper, with the exception that the new candidate set
	//is added to the field new_nbhs instead of directly replacing the out_nbh of p
	void robustPrune(tvec_point* p, parlay::sequence<pid> candidates, parlay::sequence<tvec_point*> &v, double alpha) {
    // Make sure the candidate set does not include p (done later) add
    // out neighbors of p to the candidate set.
    // std::cout << "here1" << std::endl; 
		for (size_t i=0; i<size_of(p->out_nbh); i++) {
			candidates.push_back(std::make_pair(p->out_nbh[i],
				distance(v[p->out_nbh[i]]->coordinates.begin(), p->coordinates.begin(), d)));
		}

		// Sort the candidate set in reverse order according to distance from p.
		auto less = [&] (pid a, pid b) {return a.second < b.second;};
		parlay::sort_inplace(candidates, less);

		fine_sequence new_nbhs = fine_sequence();
    // new_nbhs.reserve(maxDeg);

    size_t candidate_idx = 0;
		while (new_nbhs.size() <= maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      int p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p->id || p_star == -1) {
        continue;
      }

      new_nbhs.push_back(p_star);

      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        int p_prime = candidates[i].first;
        if (p_prime != -1) {
          float dist_starprime = distance(v[p_star]->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = -1;
          }
        }
      }
		}
		// std::cout << "here2" << std::endl;
		// std::cout << new_nbhs.size() << std::endl;  
		add_new_nbh(new_nbhs, p);
		// std::cout << "here3" << std::endl; 
	}


	//robustPrune routine as found in DiskANN paper, with the exception that the new candidate set
	//is added to the field new_nbhs instead of directly replacing the out_nbh of p
	void robustPrune(Tvec_point<T>* p, parlay::sequence<int> candidates, parlay::sequence<Tvec_point<T>*> &v, double alpha){
    parlay::sequence<pid> cc;
    cc.reserve(candidates.size() + size_of(p->out_nbh));
    for (size_t i=0; i<candidates.size(); ++i) {
      cc.push_back(std::make_pair(candidates[i], distance(v[candidates[i]]->coordinates.begin(), p->coordinates.begin(), d)));
    }
    return robustPrune(p, std::move(cc), v, alpha);
	}

	void build_index(parlay::sequence<Tvec_point<T>*> &v, bool from_empty = true, bool two_pass = true){
		clear(v);
		//populate with random edges
		// if(not from_empty) random_index(v, maxDeg);
		//find the medoid, which each beamSearch will begin from
		find_approx_medoid(v);
		build_index_inner(v, 1.0, 2, .1);
		if(two_pass) build_index_inner(v, r2_alpha, 2, .1);
	}

	void build_index_inner(parlay::sequence<Tvec_point<T>*> &v, double alpha = 1.0, double base = 2,
	double max_fraction = 1.0, bool random_order = false){
		size_t n = v.size();
		size_t inc = 0;
		parlay::sequence<int> rperm;
		if(random_order) rperm = parlay::random_permutation<int>(static_cast<int>(n), time(NULL));
		else rperm = parlay::tabulate(v.size(), [&] (int i) {return i;});
		size_t count = 0;
		while(count < n){
			size_t floor;
			size_t ceiling;
			if(pow(base,inc) <= 1000000){
				floor = static_cast<size_t>(pow(base, inc))-1;
				ceiling = std::min(static_cast<size_t>(pow(base, inc+1)), n)-1;
				count = std::min(static_cast<size_t>(pow(base, inc+1)), n)-1;
			} else{
				floor = count;
				ceiling = std::min(count + static_cast<size_t>(1000000), n)-1;
				count += static_cast<size_t>(1000000);
			}
			//search for each node starting from the medoid, then call
			//robustPrune with the visited list as its candidate set
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				parlay::sequence<pid> visited = (beam_search(v[rperm[i]], v, medoid, beamSize, d)).second;
				if(report_stats) v[rperm[i]]->cnt = visited.size();
				robustPrune(v[rperm[i]], visited, v, alpha);
			});
			//make each edge bidirectional by first adding each new edge
			//(i,j) to a sequence, then semisorting the sequence by key values
			parlay::sequence<parlay::sequence<index_pair>> to_flatten =
				parlay::sequence<parlay::sequence<index_pair>>(ceiling-floor);
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				parlay::sequence<index_pair> edges = parlay::sequence<index_pair>(size_of(v[rperm[i]]->new_nbh));
				for(int j=0; j<size_of(v[rperm[i]]->new_nbh); j++){
					edges[j] = std::make_pair((v[rperm[i]]->new_nbh)[j], rperm[i]);  
				}
				to_flatten[i-floor] = edges;
				synchronize(v[rperm[i]]);
			});
			// std::cout << "here3" << std::endl; 
			auto new_edges = parlay::flatten(to_flatten);
			auto grouped_by = parlay::group_by_key(new_edges);
			//finally, add the bidirectional edges; if they do not make
			//the vertex exceed the degree bound, just add them to out_nbhs;
			//otherwise, use robustPrune on the vertex with user-specified alpha
			// std::cout << grouped_by.size() << std::endl; 
			parlay::parallel_for(0, grouped_by.size(), [&] (size_t j){
				size_t index = grouped_by[j].first;
				parlay::sequence<int> candidates = grouped_by[j].second;
				int newsize = candidates.size() + size_of(v[index]->out_nbh);
				if(newsize <= maxDeg){
					for(const int& k : candidates) add_nbh(k, v[index]);
				} else{
					robustPrune(v[index], candidates, v, alpha);
					synchronize(v[index]);
				}
			});
			inc += 1;
		}
	}

	void searchNeighbors(parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<Tvec_point<T>*> &v, int beamSizeQ, int k){
		searchFromSingle(q, v, beamSizeQ, k, d, medoid);
	}
};