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
#include "common/geometry.h"
#include <random>
#include <set>
#include <math.h>

extern bool report_stats;


template<typename T>
struct knn_index{
	int maxDeg;
	int beamSize;
	int k; 
	double r2_alpha; //alpha parameter for round 2 of robustPrune
	unsigned d;
	using tvec_point = Tvec_point<T>;
	using fvec_point = Tvec_point<float>;
	tvec_point* medoid; 
	using pid = std::pair<int, float>;
	using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
	using index_pair = std::pair<int, int>;
	using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));

	knn_index(int md, int bs, int kk, double a, unsigned dim) : maxDeg(md), beamSize(bs), k(kk), r2_alpha(a), d(dim) {}

	void clear(parlay::sequence<tvec_point*> &v){
		size_t n = v.size();
		parlay::parallel_for(0, n, [&] (size_t i){
			parlay::sequence<int> clear_seq = parlay::sequence<int>();
			v[i]->out_nbh = clear_seq;
		});
	}

	//give each vertex maxDeg random out neighbors
	//note that this function as of now may have some contention when running in parallel
	//TODO rewrite using parlay::random functions
	void random_index(parlay::sequence<tvec_point*> &v){
		size_t n = v.size(); 
		parlay::parallel_for(0, n, [&] (size_t i){
	    	std::random_device rd;    
			std::mt19937 rng(rd());   
			std::uniform_int_distribution<int> uni(0,n-1); 

			//use a set to make sure each out neighbor is unique
			std::set<int> indexset;
			while(indexset.size() < maxDeg){
				int j = uni(rng);
				//disallow self edges
				if(j != i) indexset.insert(j);
			}
			
			for (std::set<int>::iterator it=indexset.begin(); it!=indexset.end(); ++it){
        		v[i] -> out_nbh.push_back(*it);
			} 

	    }, 1
	    );
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
				float result = c1[i] + c2[i];
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
	void robustPrune(tvec_point* p, parlay::sequence<pid> candidates, parlay::sequence<tvec_point*> &v, double alpha){
		//make sure the candidate set does not include p
		auto pred = [&] (pid a){return a.first == p->id;};
		if(find_if(candidates, pred) != candidates.end()) candidates.erase(find_if(candidates, pred));
		//add out neighbors of p to the candidate set
		for(int i=0; i<(p->out_nbh.size()); i++){
			candidates.push_back(std::make_pair(p->out_nbh[i], 
				distance(v[p->out_nbh[i]]->coordinates.begin(), p->coordinates.begin(), d)));
		}  
		//sort the candidate set in reverse order according to distance from p
		auto less = [&] (pid a, pid b) {return a.second > b.second;};
		parlay::sort_inplace(candidates, less);
		parlay::sequence<int> new_nbhs = parlay::sequence<int>();
		while(new_nbhs.size() <= maxDeg && candidates.size() > 0){
			int c = candidates.size();
			int p_star = candidates[c-1].first;
			candidates.pop_back();
			new_nbhs.push_back(p_star);
			parlay::sequence<int> to_delete = parlay::sequence<int>();
			for(int i=0; i<c-1; i++){
				int p_prime = candidates[i].first;
				float dist_starprime = distance(v[p_star]->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
				float dist_pprime = distance(p->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
				if(alpha*dist_starprime <= dist_pprime){
					to_delete.push_back(i);
				}
			}
			if(to_delete.size() > 0){
				for(int i=0; i<to_delete.size(); i++){
					candidates.erase(candidates.begin()+to_delete[i]-i);
				}
			}
		}
		p->new_out_nbh = new_nbhs;
	}

	//robustPrune routine as found in DiskANN paper, with the exception that the new candidate set
	//is added to the field new_nbhs instead of directly replacing the out_nbh of p
	void robustPrune(Tvec_point<T>* p, parlay::sequence<int> candidates, parlay::sequence<Tvec_point<T>*> &v, double alpha){
		//make sure the candidate set does not include p
		if(parlay::find(candidates, p->id) != candidates.end()) candidates.erase(parlay::find(candidates, p->id));
		//add out neighbors of p to the candidate set
		for(int i=0; i<(p->out_nbh.size()); i++){
			candidates.push_back(p->out_nbh[i]);
		}  
		//sort the candidate set in reverse order according to distance from p
		auto less = [&] (int a, int b){
			return distance(v[a]->coordinates.begin(), p->coordinates.begin(), d) > 
			distance(v[b]->coordinates.begin(), p->coordinates.begin(), d);
		};
		parlay::sort_inplace(candidates, less);
		parlay::sequence<int> new_nbhs = parlay::sequence<int>();
		while(new_nbhs.size() <= maxDeg && candidates.size() > 0){
			int c = candidates.size();
			int p_star = candidates[c-1];
			candidates.pop_back();
			new_nbhs.push_back(p_star);
			parlay::sequence<int> to_delete = parlay::sequence<int>();
			for(int i=0; i<c-1; i++){
				int p_prime = candidates[i];
				float dist_starprime = distance(v[p_star]->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
				float dist_pprime = distance(p->coordinates.begin(), v[p_prime]->coordinates.begin(), d);
				if(alpha*dist_starprime <= dist_pprime){
					to_delete.push_back(i);
				}
			}
			if(to_delete.size() > 0){
				for(int i=0; i<to_delete.size(); i++){
					candidates.erase(candidates.begin()+to_delete[i]-i);
				}
			}
		}
		p->new_out_nbh = new_nbhs;
	}

	void build_index(parlay::sequence<Tvec_point<T>*> &v, bool from_empty = true, bool random_order = true){
		clear(v);
		//populate with random edges
		if(not from_empty) random_index(v);
		//find the medoid, which each beamSearch will begin from
		find_approx_medoid(v);
		size_t n = v.size();
		size_t inc = 0;
		parlay::sequence<int> rperm;
		if(random_order) rperm = parlay::random_permutation<int>(static_cast<int>(n), time(NULL));
		else rperm = parlay::tabulate(v.size(), [&] (int i) {return i;});
		while(pow(2, inc) < n){
			size_t floor = static_cast<size_t>(pow(2, inc))-1;
			size_t ceiling = std::min(static_cast<size_t>(pow(2, inc+1)), n)-1;
			//search for each node starting from the medoid, then call
			//robustPrune with the visited list as its candidate set
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				parlay::sequence<pid> visited = (beam_search(v[rperm[i]], v, medoid, beamSize, d)).second;
				if(report_stats) v[rperm[i]]->cnt = visited.size();
				robustPrune(v[rperm[i]], visited, v, 1);
			});
			//make each edge bidirectional by first adding each new edge
			//(i,j) to a sequence, then semisorting the sequence by key values
			parlay::sequence<parlay::sequence<index_pair>> to_flatten = 
				parlay::sequence<parlay::sequence<index_pair>>(ceiling-floor);
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				parlay::sequence<int> new_nbh = v[rperm[i]]->new_out_nbh;
				parlay::sequence<index_pair> edges = parlay::sequence<index_pair>(new_nbh.size());
				for(int j=0; j<new_nbh.size(); j++){
					edges[j] = std::make_pair(new_nbh[j], rperm[i]);
				}
				to_flatten[i-floor] = edges;
				v[rperm[i]]->out_nbh = new_nbh;
				(v[rperm[i]]->new_out_nbh).clear();
			});
			auto new_edges = parlay::flatten(to_flatten);
			auto grouped_by = parlay::group_by_key(new_edges);
			//finally, add the bidirectional edges; if they do not make 
			//the vertex exceed the degree bound, just add them to out_nbhs;
			//otherwise, use robustPrune on the vertex with user-specified alpha
			parlay::parallel_for(0, grouped_by.size(), [&] (size_t j){
				size_t index = grouped_by[j].first;
				parlay::sequence<int> candidates = grouped_by[j].second;
				int newsize = candidates.size() + (v[index]->out_nbh).size();
				if(newsize <= maxDeg){
					for(int k=0; k<candidates.size(); k++){ //try concatenating instead of pushing back
						(v[index]->out_nbh).push_back(candidates[k]);
					}
				} else{ 
					robustPrune(v[index], candidates, v, r2_alpha);
					parlay::sequence<int> new_nbh = v[index]->new_out_nbh;
					v[index]->out_nbh = new_nbh;
					(v[index]->new_out_nbh).clear();
				}
			});
			inc += 1; 
		}
	}

	void searchNeighbors(parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<Tvec_point<T>*> &v){
		if((k+1)>beamSize){
			std::cout << "Error: beam search parameter L = " << beamSize << " same size or smaller than k = " << k << std::endl;
			abort();
		}
		parlay::parallel_for(0, q.size(), [&] (size_t i){
			parlay::sequence<int> neighbors = parlay::sequence<int>(k);
			parlay::sequence<pid> beamElts = (beam_search(q[i], v, medoid, beamSize, d)).first;
			//the first element of the frontier may be the point itself
			//if this occurs, do not report it as a neighbor
			if(beamElts[0].first==i){
				for(int j=0; j<k; j++){
					neighbors[j] = beamElts[j+1].first;
				}
			} else{
				for(int j=0; j<k; j++){
					neighbors[j] = beamElts[j].first;
				}
			}
			q[i]->ngh = neighbors;
		});
	}
};