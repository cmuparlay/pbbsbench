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
#include <set>

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

template<class fvec_point>
std::pair<parlay::sequence<int>, std::set<int>> beam_search(fvec_point* p, parlay::sequence<fvec_point*> &v, 
	fvec_point* medoid, int beamSize){
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