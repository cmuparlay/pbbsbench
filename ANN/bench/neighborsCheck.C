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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "benchUtils.h"

using namespace benchIO;

//start with 1@1, later move to other types
int checkNeighbors(int k, parlay::sequence<fvec_point> groundTruth, parlay::sequence<long> neighbors){
	int numCorrect = 0;
	size_t n = (neighbors.size())/(k+1);
	for(int i=0; i<n; i++){
		int reported_index = neighbors[(k+1)*i+1];
		int true_index = (groundTruth[i].coordinates)[0];
		if(reported_index == true_index) numCorrect += 1;
	}
	float recall = static_cast<float>(numCorrect)/static_cast<float>(n);
	std:: cout << "Recall 1@1: " << recall << std::endl; 
	return 0;
}

int main(int argc, char* argv[]) {
	commandLine P(argc, argv, "[-k {1,...,100}] <inFile> <outfile>");
	pair<char*,char*> fnames = P.IOFileNames();
	char* iFile = fnames.first; //the ground truth
	char* oFile = fnames.second; //the output of the algorithm

	int k = P.getOptionIntValue("-k",1);
	if (k > 100 || k < 1) P.badArgument();

	parlay::sequence<long> neighbors = readIntSeqFromFile<long>(oFile);
	auto groundTruth = parse_fvecs(iFile);
	checkNeighbors(k, groundTruth, neighbors);
	std::cout << "Executing check file" << std::endl; 

}
