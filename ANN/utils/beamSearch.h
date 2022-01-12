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

#ifndef BEAMSEARCH
#define BEAMSEARCH

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include <set>

using pid = std::pair<int, float>;

//returns true if F \setminus V = emptyset 
bool intersect_nonempty(parlay::sequence<pid> &V, parlay::sequence<pid> &F){
	for(int i=0; i<F.size(); i++){
		auto pred = [&] (pid a) {return F[i].first == a.first;};
		if(parlay::find_if(V, pred) == V.end()) return true; 
	} 
	return false;
}

//will only be used when there is an element in F that is not in V
//hence the ``return 0" line will never be called
pid id_next(parlay::sequence<pid> &V, parlay::sequence<pid> &F){
	for(int i=0; i<F.size(); i++){
		auto pred = [&] (pid a) {return F[i].first == a.first;};
		if(parlay::find_if(V, pred) == V.end()) return F[i];
	}
	return std::make_pair(0,0);
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

//takes in two sorted sequences and returns a sorted union
template<class fvec_point>
parlay::sequence<pid> seq_union(parlay::sequence<pid> &P, parlay::sequence<pid> &Q, fvec_point* p, 
	parlay::sequence<fvec_point*> &v, unsigned d){
	auto less = [&] (pid a, pid b){
		return a.second < b.second;
	};
	pid* first1 = P.begin();
	pid* last1 = P.end();
	pid* first2 = Q.begin();
	pid* last2 = Q.end();
	parlay::sequence<pid> result = parlay::sequence<pid>();
	while(true){
		if(first1 == last1){
			while(first2 != last2){
				result.push_back(*first2);
				++first2;
			}
			return result;
		} else if(first2 == last2){
			while(first1 != last1){
				result.push_back(*first1);
				++first1;
			}
			return result;
		} 
		if(less(*first1, *first2)){
			result.push_back(*first1);
			++first1;
		} else if(less(*first2, *first1)){
			result.push_back(*first2);
			++first2;
		} else{
			if(first1->first == first2->first){
				result.push_back(*first1);
				++first1;
				++first2; 
			} else{
				result.push_back(*first1);
				result.push_back(*first2);
				++first1;
				++first2;
			}
		}
	}
	return result;
}

template<class fvec_point>
std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(fvec_point* p, parlay::sequence<fvec_point*> &v, 
	fvec_point* medoid, int beamSize, unsigned d){
	//initialize data structures
	parlay::sequence<pid> visited = parlay::sequence<pid>();
	parlay::sequence<pid> frontier = parlay::sequence<pid>();
	auto less = [&] (pid a, pid b){
		return a.second < b.second;
	};
	//the frontier starts with the medoid
	frontier.push_back(std::make_pair(medoid->id, distance(medoid, p, d)));
	//terminate beam search when the entire frontier has been visited
	while(intersect_nonempty(visited, frontier)){
		//the next node to visit is the unvisited frontier node that is closest to p
		pid currentPid = (id_next(visited, frontier));
		fvec_point* current = v[currentPid.first]; 
		auto g = [&] (int a){
			auto pred = [&] (pid p) {return p.first == a;};
			if(parlay::find_if(visited, pred) == visited.end()) return true;
			return false;   
		};
		auto candidates = parlay::filter(current->out_nbh, g);
		parlay::sequence<pid> pairCandidates = parlay::sequence<pid>(candidates.size());
		for(int i=0; i<candidates.size(); i++) pairCandidates[i] = std::make_pair(candidates[i], distance(v[candidates[i]], p, d));
		auto sortedCandidates = parlay::sort(pairCandidates, less);
		frontier = seq_union(frontier, sortedCandidates, p, v, d);
		if(frontier.size() > beamSize) frontier.erase(frontier.begin()+beamSize, frontier.end());
		//add the node to the visited list
		visited.push_back(currentPid);
	} 
	// parlay::sequence<int> frontier = parlay::sequence<int>(pairFrontier.size());
	// for(int i=0; i<pairFrontier.size(); i++) frontier[i] = (pairFrontier[i]).first;
	return std::make_pair(frontier, visited);
}

#endif