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
#include <math.h>


template<class fvec_point>
struct knn_index{
	int maxDeg;
	int beamSize;
	int k; 
	double r2_alpha; //alpha parameter for round 2 of robustPrune
	fvec_point* medoid; 
	using slice_fvec = decltype(make_slice(parlay::sequence<fvec_point*>()));
	using index_pair = std::pair<int, int>;
	using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));

	knn_index(int md, int bs, int kk, double a) : maxDeg(md), beamSize(bs), k(kk), r2_alpha(a) {}

	void clear(parlay::sequence<fvec_point*> &v){
		size_t n = v.size();
		parlay::parallel_for(0, n, [&] (size_t i){
			parlay::sequence<int> clear_seq = parlay::sequence<int>();
			v[i]->out_nbh = clear_seq;
		});
	}

	//give each vertex maxDeg random out neighbors
	void random_index(parlay::sequence<fvec_point*> &v){
		size_t n = v.size(); 
		parlay::parallel_for(0, n, [&] (size_t i){
	    	std::random_device rd;    
			std::mt19937 rng(rd());   
			std::uniform_int_distribution<int> uni(0,n-1); 

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
	void find_approx_medoid(parlay::sequence<fvec_point*> &v){
		size_t n = v.size();
		int d = ((v[0]->coordinates).size())/4;
		parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v), d);
		fvec_point centroidp = typename fvec_point::fvec_point();
		centroidp.coordinates = parlay::make_slice(centroid);
		medoid = medoid_helper(&centroidp, parlay::make_slice(v), d); 
	}

	//returns true if F \setminus V = emptyset 
	bool intersect_nonempty(std::set<int> V, parlay::sequence<int> F){
		int Fsize = F.size();
		for(int i=0; i<Fsize; i++){
			if(V.find(F[i]) == V.end()) return true; 
		}
		return false;
	}

	//will only be used when there is an element in F that is not in V
	//hence the ``return 0" line will never be called
	int id_next(std::set<int> V, parlay::sequence<int> F){
		int fsize = F.size();
		for(int i=0; i<fsize; i++){
			if(V.find(F[i]) == V.end()) return F[i];
		}
		return 0;
	}

	//for debugging
	void print_seq(parlay::sequence<int> seq){
		int fsize = seq.size();
		std::cout << "["; 
		for(int i=0; i<fsize; i++){
			std::cout << seq[i] << ", ";
		}
		std::cout << "]" << std::endl; 
	}

	void print_set(std::set<int> myset){
		std::cout << "["; 
		for (std::set<int>::iterator it=myset.begin(); it!=myset.end(); ++it){
			std::cout << *it << ", ";
		}
		std::cout << "]" << std::endl; 
	}


	std::pair<parlay::sequence<int>, std::set<int>> beam_search(fvec_point* p, parlay::sequence<fvec_point*> &v){
		//initialize data structures
		std::set<int> visited;
		std::set<int> frontier; 
		parlay::sequence<int> sortedFrontier = parlay::sequence<int>();
		//the frontier starts with the medoid
		sortedFrontier.push_back(medoid->id);
		frontier.insert(medoid->id);
		//terminate beam search when the entire frontier has been visited
		while(intersect_nonempty(visited, sortedFrontier)){
			//the next node to visit is the unvisited frontier node that is closest to p
			int currentIndex = id_next(visited, sortedFrontier);
			fvec_point* current = v[currentIndex]; 
			parlay::sequence<int> outnbh = current->out_nbh;
			int outsize = outnbh.size();
			//add the outneighbors of the visited node to the frontier if they are not already in it
			for(int i=0; i<(outsize); i++){
				if(frontier.find(outnbh[i]) == frontier.end()){
					auto less = [&] (int a){
						return distance(v[a], p) < distance(v[outnbh[i]], p);
					};
					int insertion_point = parlay::internal::binary_search(parlay::make_slice(sortedFrontier), less);
					const int to_insert = outnbh[i];
					sortedFrontier.insert(sortedFrontier.begin()+insertion_point, to_insert);
					frontier.insert(outnbh[i]); 
				}
			}
			//remove nodes from frontier if too large
			if(sortedFrontier.size() > beamSize){ //beamSize is a global var taken from command line
				for(int i=sortedFrontier.size(); i>beamSize; i--){
					int to_del = sortedFrontier[i-1];
					sortedFrontier.pop_back();
					frontier.erase(to_del);
				}
			}
			//add the node to the visited list
			visited.insert(current->id);
		} 
		return std::make_pair(sortedFrontier, visited);
	}

	//robustPrune routine as found in DiskANN paper, with the exception that the new candidate set
	//is added to the field new_nbhs instead of directly replacing the out_nbh of p
	void robustPrune(fvec_point* p, std::set<int> candidates, parlay::sequence<fvec_point*> &v, double alpha){
		//add out neighbors of p to the candidate set
		for(int i=0; i<(p->out_nbh.size()); i++){
			candidates.insert(p->out_nbh[i]);
		}
		//transform to a sequence so that it can be sorted
		parlay::sequence<int> Candidates = parlay::sequence<int>();
		for (std::set<int>::iterator it=candidates.begin(); it!=candidates.end(); ++it){
        	Candidates.push_back(*it);
		} 
		//sort the candidate set in reverse order according to distance from p
		auto less = [&] (int a, int b){
			return distance(v[a], p) > distance(v[b], p);
		};
		auto sortedCandidates = parlay::sort(Candidates, less);
		parlay::sequence<int> new_nbhs = parlay::sequence<int>();
		while(new_nbhs.size() <= maxDeg && sortedCandidates.size() > 0){
			int c = sortedCandidates.size();
			int p_star = sortedCandidates[c-1];
			sortedCandidates.pop_back();
			new_nbhs.push_back(p_star);
			parlay::sequence<int> to_delete = parlay::sequence<int>();
			for(int i=0; i<c-1; i++){
				int p_prime = sortedCandidates[i];
				float dist_starprime = distance(v[p_star], v[p_prime]);
				float dist_pprime = distance(p, v[p_prime]);
				if(alpha*dist_starprime <= dist_pprime){
					to_delete.push_back(i);
				}
			}
			if(to_delete.size() > 0){
				for(int i=0; i<to_delete.size(); i++){
					sortedCandidates.erase(sortedCandidates.begin()+to_delete[i]-i);
				}
			}
		}
		p->new_out_nbh = new_nbhs;
	}

	void build_index(parlay::sequence<fvec_point*> &v){
		clear(v);
		random_index(v);
		find_approx_medoid(v);
		size_t n = v.size();
		size_t inc = 0;
		while(pow(2, inc) < n){
			size_t floor = static_cast<size_t>(pow(2, inc))-1;
			size_t ceiling = std::min(static_cast<size_t>(pow(2, inc+1)), n)-1;
			//search for each node starting from the medoid, then call
			//robustPrune with the visited list as its candidate set
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				robustPrune(v[i], (beam_search(v[i], v)).second, v, 1);
			});
			//make each edge bidirectional by first adding each new edge
			//(i,j) to a sequence, then semisorting the sequence by key values
			parlay::sequence<parlay::sequence<index_pair>> to_flatten = 
				parlay::sequence<parlay::sequence<index_pair>>(ceiling-floor);
			parlay::parallel_for(floor, ceiling, [&] (size_t i){
				parlay::sequence<int> new_nbh = v[i]->new_out_nbh;
				parlay::sequence<index_pair> edges = parlay::sequence<index_pair>(new_nbh.size());
				for(int j=0; j<new_nbh.size(); j++){
					edges[j] = std::make_pair(new_nbh[j], i);
				}
				to_flatten[i-floor] = edges;
				v[i]->out_nbh = new_nbh;
				v[i]->new_out_nbh = parlay::sequence<int>();
			});
			auto new_edges = parlay::flatten(to_flatten);
			auto grouped_by = parlay::group_by_key(new_edges);
			//finally, add the bidirectional edges; if they do not make 
			//the vertex exceed the degree bound, just add them to out_nbhs;
			//otherwise, use robustPrune on the vertex with user-specified alpha
			parlay::parallel_for(0, grouped_by.size(), [&] (size_t i){
				size_t index = grouped_by[i].first;
				parlay::sequence<int> candidates = grouped_by[i].second;
				int newsize = candidates.size() + (v[index]->out_nbh).size();
				if(newsize <= maxDeg){
					for(int i=0; i<candidates.size(); i++){
						(v[index]->out_nbh).push_back(candidates[i]);
					}
				} else{
					std::set<int> candidateSet;
					for(int i=0; i<candidates.size(); i++){
						candidateSet.insert(candidates[i]);
					}
					robustPrune(v[index], candidateSet, v, r2_alpha);
					parlay::sequence<int> new_nbh = v[index]->new_out_nbh;
					v[index]->out_nbh = new_nbh;
					v[index]->new_out_nbh = parlay::sequence<int>();
				}
			});
			inc += 1; 
		}
	}
};