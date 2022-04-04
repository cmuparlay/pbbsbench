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


template<typename T>
void clear(parlay::sequence<Tvec_point<T>*> &v){
	size_t n = v.size();
	parlay::parallel_for(0, n, [&] (size_t i){
		v[i]->out_nbh.clear();
	});
}   

//TODO fix to use custom allocator
template<typename T>
void random_index(parlay::sequence<Tvec_point<T>*> &v, int maxDeg){
	size_t n = v.size(); 
	parlay::parallel_for(0, n, [&] (size_t i){
		// std::cout << "here1" << std::endl; 
    	std::random_device rd;    
		std::mt19937 rng(rd());   
		std::uniform_int_distribution<int> uni(0,n-1); 
    	// std::cout << "here2" << std::endl; 
		//use a set to make sure each out neighbor is unique
		std::set<int> indexset;
		while(indexset.size() < maxDeg){
			int j = uni(rng);;
			//disallow self edges
			if(j != i) indexset.insert(j);
		}
		// std::cout << "here3" << std::endl; 
		for (std::set<int>::iterator it=indexset.begin(); it!=indexset.end(); ++it){
    		v[i] -> out_nbh.push_back(*it);
		} 

    }, 1
    );
}