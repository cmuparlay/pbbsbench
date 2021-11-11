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


template<class fvec_point>
struct knn_index{
	int maxDeg;
	int beamSize;
	int k; 
	fvec_point* medoid; 
	using slice_fvec = decltype(make_slice(parlay::sequence<fvec_point*>()));

	knn_index(int md, int bs, int kk) : maxDeg(md), beamSize(bs), k(kk) {}

	//give each vertex maxDeg random out neighbors
	void random_index(parlay::sequence<fvec_point*> v){
		size_t n = v.size(); 
		parlay::parallel_for(0, n, [&] (size_t i){
	    	std::random_device rd;    
			std::mt19937 rng(rd());   
			std::uniform_int_distribution<int> uni(0,n); 

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
	void find_approx_medoid(parlay::sequence<fvec_point*> v){
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

	int id_next(std::set<int> V, parlay::sequence<int> F){
		int fsize = F.size();
		for(int i=0; i<fsize; i++){
			if(V.find(F[i]) == V.end()) return F[i];
		}
		return 0;
	}


	//actually should redo this using sets, just also need a sorted list for the frontier
	//so frontier has a sorted list AND a set, while visited has just a set
	std::pair<std::set<int>, std::set<int>> beam_search(fvec_point* p, parlay::sequence<fvec_point*> v){
		std::set<int> visited;
		std::set<int> frontier; 
		parlay::sequence<int> sortedFrontier = parlay::sequence<int>();
		sortedFrontier.push_back(medoid->id);
		frontier.insert(medoid->id);
		while(intersect_nonempty(visited, sortedFrontier)){
			//identify the node we are visiting
			int currentIndex = id_next(visited, sortedFrontier);
			fvec_point* current = v[currentIndex]; 
			parlay::sequence<int> outnbh = current->out_nbh;
			for(int i=0; i++; i<(outnbh.size())){
				if(frontier.find(outnbh[i]) != frontier.end()){
					auto less = [&] (int a){
						return distance(v[a], p) < distance(v[outnbh[i]], p);
						// return a<outnbh[i]; 
					};
					int insertion_point = parlay::internal::binary_search(parlay::make_slice(sortedFrontier), less);
					// sortedFrontier.begin()+insertion_point; 
					const int to_insert = outnbh[i];
					sortedFrontier.insert(sortedFrontier.begin()+insertion_point, to_insert);
					frontier.insert(outnbh[i]);
				}
			}
			//remove nodes from frontier if too large
			if(sortedFrontier.size() > beamSize){ //beamSize is a global var taken from command line
				for(int i=sortedFrontier.size(); i>beamSize; i--){
					int to_del = sortedFrontier[i];
					sortedFrontier.pop_back();
					frontier.erase(to_del);
				}
			}
			//add the node to the visited list
			visited.insert(current->id);
		}
		return std::make_pair(frontier, visited);
	}

	// void robustPrune(fvec_point* p, std::set<int> candidates){
	// 	//transform the candidate set to a sequence, add the out_nbhs of p, and sort
	// 	for(int i=0; i<(p->out_nbh.size()); i++){
	// 		candidates.insert(p->out_nbh[i]);
	// 	}
	// 	parlay::sequence<int> sortedCandidates = parlay::sequence<int>(candidates.size());
	// 	for(int i=0; i<candidates.size(); i++){ //need to use the set iterator

	// 	}
	// }

	void build_index(parlay::sequence<fvec_point*> v){
		random_index(v);
		find_approx_medoid(v);
		beam_search(v[0], v);
	}
};