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

#include <iostream>
#include "float.h"
#include <algorithm>
#include <cstring>
#include <set>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "benchUtils.h"


using namespace benchIO;

//make a set of the real neighbors
//cut seq of the reported neighbors
//loop over seq, add one if find() reports something
int checkNeighbors(int k, parlay::sequence<ivec_point> groundTruth, parlay::sequence<long> neighbors, parlay::sequence<int> recall_types){
	size_t n = (neighbors.size())/(k+1); 
	for(int j=0; j<recall_types.size(); j++){
		int r = recall_types[j];
		int numCorrect = 0;
		for(int i=0; i<n; i++){
			std::set<int> reported_nbhs;
			for(int l=0; l<r; l++) reported_nbhs.insert(neighbors[i*(k+1)+1+l]);
			for(int l=0; l<r; l++)
			  if (reported_nbhs.find((groundTruth[i].coordinates)[l]) != reported_nbhs.end())
			    numCorrect += 1;
		}
		float recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
		std:: cout << "Recall " << r << "@" << r << ": " << recall << std::endl; 
	}
	return 1;       
}

//parses a string of the form "[i, j, l]" to a vector of ints
//removes all ints i which are greater than k
parlay::sequence<int> parse_recall(std::string recall, int k){
	parlay::sequence<int> recall_vec = parlay::sequence<int>();
	for(std::string::size_type i = 0; i < recall.size(); i++){
		char c = recall[i];
			if(c == ' ' | c == '[' | c == ']' | c == ',') continue;
			else{
				std::string s = std::string(1, c);
				char next = recall[i+1];
				while(next != ' ' && next != '[' && next != ']' && next != ','){
					s += next;
					i += 1;
					next = recall[i+1];
				}
				int elt = std::stoi(s);
				if(elt <= k) recall_vec.push_back(elt);
			}
	}
	return recall_vec;
}

int main(int argc, char* argv[]) {
	commandLine P(argc, argv, "[-r <recall> <inFile> <outfile>");
	pair<char*,char*> fnames = P.IOFileNames();
	char* iFile = fnames.first; //the ground truth
	char* oFile = fnames.second; //the output of the algorithm
	std::string recall = P.getOptionValue("-r", "[1]");
	parlay::sequence<long> neighbors = readIntSeqFromFile<long>(oFile);    
	auto groundTruth = parse_ivecs(iFile);
	int k = (neighbors.size())/(groundTruth.size())-1;
	parlay::sequence<int> recall_vec = parse_recall(recall, k);
	return checkNeighbors(k, groundTruth, neighbors, recall_vec);
}