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

#ifndef _myRMQ_hpp_
#define _myRMQ_hpp_

#include <iostream>
#include "parallel.h"
#include "utils.h"
#include "math.h"
#define BSIZE 128
using namespace std;

class myRMQ {
protected:
  uintT* a;
  uintT** table;
  long n, m, depth;

 public:
  myRMQ(uintT* _a, long _n){
    a = _a;
    n = _n;
    m = 1 + (n-1)/BSIZE;
    precomputeQueries();
  }

  void precomputeQueries() {
    depth = utils::log2Up(m+1);
    table = newA(uintT*, depth);
    parallel_for(long k=0; k<depth; k++)
      table[k] = newA(uintT, m);
    parallel_for(long i=0; i < m; i++) {
      long start = i*BSIZE;
      long end = min(start+BSIZE,n);
      long k = i*BSIZE;
      for (long j = start+1; j < end; j++) 
	if (a[j] < a[k]) k = j;
      table[0][i] = k;
    }

    long dist = 1;
    for (long j=1; j<depth; j++) {
      parallel_for (long i=0; i < m-dist; i++) {
	if (a[table[j-1][i]] <= a[table[j-1][i+dist]])
	  table[j][i] = table[j-1][i];
	else table[j][i] = table[j-1][i+dist];
      }
      parallel_for (long l = m-dist; l<m; l++)
	table[j][l] = table[j-1][l];
      dist*=2;
    }
  }

  long query(long i, long j) {
    //same block
    if (j-i < BSIZE) {
      long r = i;
      for (long k = i+1; k <= j; k++) 
	if (a[k] < a[r]) r = k;
      return r;
    } 
    long block_i = i/BSIZE;
    long block_j = j/BSIZE;
    long min = i;
    for (long k=i+1; k < (block_i+1)*BSIZE; k++)
      if(a[k] < a[min]) min = k;
    for (long k=j; k>=(block_j)*BSIZE; k--)
      if(a[k] < a[min]) min = k;
    if (block_j == block_i + 1) return min;
    long outOfBlockMin;

    //not same or adjacent blocks
    if (block_j > block_i + 1) {
      block_i++;
      block_j--;
      if (block_j == block_i) 
	outOfBlockMin = table[0][block_i];
      else if(block_j == block_i + 1) 
	outOfBlockMin = table[1][block_i];
      else {
	long k = utils::log2Up(block_j - block_i+1)-1;
	long p = 1<<k; //2^k
	outOfBlockMin = a[table[k][block_i]] <= a[table[k][block_j+1-p]]
	  ? table[k][block_i] : table[k][block_j+1-p];
      }
    }
    return a[min] < a[outOfBlockMin] ? min : outOfBlockMin;
  }

  ~myRMQ() {
    parallel_for(long i=0;i<depth;i++) free(table[i]);
    free(table);
  }
};

#endif
