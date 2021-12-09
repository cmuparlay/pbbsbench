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
bool intersect_nonempty(parlay::sequence<int> &V, parlay::sequence<int> &F){
	int Fsize = F.size();
	for(int i=0; i<Fsize; i++){
		if(parlay::find(V, F[i]) == V.end()) return true;
	}
	return false;
}

//will only be used when there is an element in F that is not in V
//hence the ``return 0" line will never be called
int id_next(parlay::sequence<int> &V, parlay::sequence<int> &F){
	int fsize = F.size();
	for(int i=0; i<fsize; i++){
		if(parlay::find(V, F[i]) == V.end()) return F[i];
	}
	return 0;
}

template<class fvec_point>
std::pair<parlay::sequence<int>, parlay::sequence<int>> beam_search(fvec_point* p, parlay::sequence<fvec_point*> &v,
	fvec_point* medoid, int beamSize, unsigned d){
	//initialize data structures
	parlay::sequence<int> visited = parlay::sequence<int>();
	parlay::sequence<int> frontier = parlay::sequence<int>();
	//the frontier starts with the medoid
	frontier.push_back(medoid->id);
	//terminate beam search when the entire frontier has been visited
	while(intersect_nonempty(visited, frontier)){
		//the next node to visit is the unvisited frontier node that is closest to p
		int currentIndex = id_next(visited, frontier);
		fvec_point* current = v[currentIndex];
		//add the outneighbors of the visited node to the frontier if they are not already in it
		for(int i=0; i<((current->out_nbh).size()); i++){
			if(parlay::find(frontier, (current->out_nbh)[i]) == frontier.end()){
				auto less = [&] (int a){
					return distance(v[a], p, d) < distance(v[(current->out_nbh)[i]], p, d);
				};
				int insertion_point = parlay::internal::binary_search(parlay::make_slice(frontier), less);
				const int to_insert = (current->out_nbh)[i];
				frontier.insert(frontier.begin()+insertion_point, to_insert);
			}
		}
		if(frontier.size() > beamSize) frontier.erase(frontier.begin()+beamSize, frontier.end());
		//add the node to the visited list
		visited.push_back(current->id);
	}
	return std::make_pair(frontier, visited);
}