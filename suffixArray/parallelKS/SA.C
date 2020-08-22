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

// This is a parallel version of the algorithm described in
//  Juha Karkkainen and Peter Sanders.
//  Simple linear work suffix array construction.
//  Proc. ICALP 2003.  pp 943
// It includes code for finding the LCP
//   Written by Guy Blelloch and Julian Shun

#include <iostream>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/get_time.h"
//#include "common/range_min.h"
using namespace std;

using uchar = unsigned char;
using uint = unsigned int;
using uintPair = pair<uint,uint>;
using longInt = unsigned __int128;

inline bool leq(uint a1, uint a2,   uint b1, uint b2) {
  return(a1 < b1 || (a1 == b1 && a2 <= b2)); 
}                                                  

inline bool leq(uint a1, uint a2, uint a3, uint b1, uint b2, uint b3) {
  return(a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3))); 
}

// inline long computeLCP(uint* LCP12, uint* rank, myRMQ & RMQ, 
// 		      long j, long k, uint* s, long n){

 
//   long rank_j=rank[j]-2;
//   long rank_k=rank[k]-2;
//   if(rank_j > rank_k) {swap(rank_j,rank_k);} //swap for RMQ query

//   long l = ((rank_j == rank_k-1) ? LCP12[rank_j] 
// 	   : LCP12[RMQ.query(rank_j,rank_k-1)]);

//   long lll = 3*l;
//   if (s[j+lll] == s[k+lll]) {
//     if (s[j+lll+1] == s[k+lll+1]) return lll + 2;
//     else return lll + 1;
//   } 
//   return lll;
// }

timer t("KS SA", false);

// This recursive version requires s[n]=s[n+1]=s[n+2] = 0
// K is the maximum value of any element in s
parlay::sequence<uint> suffixArrayRec(parlay::sequence<uint> const &s,
				      size_t K, bool findLCPs) {
  t.next("head");
  size_t n = s.size() - 3; // padded with 3 nulls
  n = n+1;
  size_t n0=(n+2)/3, n1=(n+1)/3, n12=n-n0;
  auto sorted12 = parlay::sequence<uint>::uninitialized(n12);
  auto name12 = parlay::sequence<uint>::uninitialized(n12);
  auto get_first = [&] (uintPair p) {return p.first;};
  size_t bits = parlay::log2_up(K);
  //cout << n << ", " << bits << endl;
  
  if (3*bits <= 8*sizeof(uint)) {
    // if 3 chars fit into a uint then just do one radix sort
    auto C = parlay::tabulate(n12, [&] (size_t i) -> uintPair {
      size_t j = 1+(3*i)/2;
      return std::make_pair((((uint) s[j]) << 2*bits) +
			    (((uint) s[j+1]) << bits) +
			    s[j+2],
			    j);
      });

    parlay::internal::integer_sort_inplace(make_slice(C), get_first, 3*bits);
    
    // copy sorted results into sorted12
    parlay::parallel_for (1, n12, [&] (size_t i) {
	sorted12[i] = C[i].second;
	name12[i] = (C[i].first != C[i-1].first);
      });
    name12[0] = 1;
    sorted12[0] = C[0].second;
    t.next("int sort");
  } else {
    // otherwise do a comparison sort on 128 bit integers
    // with three characters and index packed in
    auto C = parlay::tabulate(n12, [&] (size_t i) -> longInt {
	size_t j = 1+(3*i)/2;
	longInt r = ((((longInt) s[j]) << 2*bits) +
		     (((size_t) s[j+1]) << bits) +
		     s[j+2]);
	return (r << 32) + j;
      });

    parlay::internal::sample_sort_inplace(make_slice(C), std::less<longInt>());

    // copy sorted results into sorted12
    longInt mask = ((((longInt) 1) << 32)-1);
    parlay::parallel_for (1, n12, [&] (size_t i) {
	sorted12[i] = C[i] & mask;
	name12[i] = (C[i] >> 32) != (C[i-1] >> 32);
      });
    name12[0] = 1;
    sorted12[0] = C[0] & mask;
    t.next("comp sort");
  }

  size_t num_names = parlay::scan_inclusive_inplace(name12, parlay::addm<uint>());
  //cout << "here c: " << name12[0] << endl;
  parlay::sequence<uint> SA12;
  //parlay::sequence<uint> LCP12;
  // recurse if names are not yet unique
  if (num_names < n12) {
    auto s12 = parlay::sequence<uint>::uninitialized(n12 + 3);
    s12[n12] = s12[n12+1] = s12[n12+2] = 0; // pad with 3 nulls

    // move mod 1 suffixes to bottom half and and mod 2 suffixes to top
    parlay::parallel_for (0, n12, [&] (size_t i) {
	uint div3 = sorted12[i]/3;
	if (sorted12[i]-3*div3 == 1) s12[div3] = name12[i];
	else s12[div3+n1] = name12[i];
      });
    name12.clear();
    sorted12.clear();

    t.next("to rec");
    SA12 = suffixArrayRec(s12, num_names+1, findLCPs);
    //SA12 = SA12_LCP.first;
    //LCP12 = SA12_LCP.second;
    s12.clear();
    
    // restore proper indices into original array
    parlay::parallel_for (0, n12, [&] (size_t i) {
	size_t l = SA12[i]; 
	SA12[i] = (l<n1) ? 3*l+1 : 3*(l-n1)+2;
      });
  } else {
    name12.clear();
    SA12 = std::move(sorted12); // suffix array is sorted array
    // if (findLCPs) {
    //   LCP12 = newA(uint, n12+3);
    //   parallel_for(long i=0; i<n12+3; i++) 
    // 	LCP12[i] = 0; //LCP's are all 0 if not recursing
    // }
  }

  // place ranks for the mod12 elements in full length array
  // mod0 locations of rank will contain garbage
  t.next("to rank");
  auto rank = parlay::sequence<uint>::uninitialized(n+2);
  rank[n]=1; rank[n+1] = 0;
  parlay::parallel_for (0, n12, [&] (size_t i) {rank[SA12[i]] = i+2;});
  t.next("set rank");
  // stably sort the mod 0 suffixes 
  // uses the fact that we already have the tails sorted in SA12
  auto mod3is1 = [&] (size_t i) {return i%3 == 1;};
  parlay::sequence<uint> s0 = parlay::filter(SA12, mod3is1);
  size_t s0n = s0.size();
  t.next("filter");
  auto D = parlay::sequence<uintPair>::uninitialized(n0);
  D[0] = make_pair(s[n-1], n-1);
  parlay::parallel_for (0, s0n, [&] (size_t i) {
      D[i+n0-s0n] = make_pair(s[s0[i]-1], s0[i]-1);});
  s0.clear();
  t.next("to sort");
  parlay::internal::integer_sort_inplace(make_slice(D), get_first, bits);
  t.next("sort");
  auto SA0 = parlay::tabulate(n0, [&] (size_t i) -> uint {return D[i].second;});
  D.clear();

  auto less = [&] (size_t i, size_t j) {
    if (i%3 == 1 || j%3 == 1) 
      return (make_pair(s[i], rank[i+1]) <
	      make_pair(s[j], rank[j+1]));
    else
      return (make_tuple(s[i], s[i+1], rank[i+2]) <
	      make_tuple(s[j], s[j+1], rank[j+2]));
  };

  //parlay::sequence<uint> SA(n-1);
  int offset = (n%3 == 1);
  parlay::sequence<uint> SA = parlay::merge(SA0.cut(offset, n0),
					  SA12.cut(1 - offset, n12),
					  less);

  t.next("merge");
  //parlay::sequence<uint> LCP;

  //get LCP from LCP12
  // if(findLCPs){
  //   LCP = newA(uint, n);  
  //   LCP[n-1] = LCP[n-2] = 0; 
  //   LCPtime.start();
  //   myRMQ RMQ(LCP12, n12+3); //simple rmq
  //   parallel_for(long i=0;i<n-2;i++){ 
  //     long j = SA[i];
  //     long k = SA[i+1];
  //     uint CLEN = 16;
  //     long ii;
  //     for (ii=0; ii < CLEN; ii++) 
  // 	if (s[j+ii] != s[k+ii]) break;
  //     if (ii != CLEN) LCP[i] = ii;
  //     else {
  //     	if (j%3 != 0 && k%3 != 0)  
  // 	  LCP[i] = computeLCP(LCP12, rank, RMQ, j, k, s, n); 
  // 	else if (j%3 != 2 && k%3 != 2)
  // 	  LCP[i] = 1 + computeLCP(LCP12, rank, RMQ, j+1, k+1, s, n);
  // 	else 
  // 	  LCP[i] = 2 + computeLCP(LCP12, rank, RMQ, j+2, k+2, s, n);
  // 	  }
  //   }
  //   LCPtime.stop();
  //   free(LCP12);
  // }
  return SA;
}

parlay::sequence<unsigned int> suffixArray(parlay::sequence<unsigned char> const &s) {
  //cout << "here 0" << endl;
  bool findLCPs = false;
  size_t n = s.size();
  auto ss = parlay::tabulate(n + 3, [&] (size_t i) -> uint {
      if (i >= n) return 0; // pad with three nulls
      else return ((uint) s[i]) + (uint) 1;
    });
  size_t k = 1 + parlay::reduce(ss, parlay::maxm<uint>());

  //cout << "here" << endl;
  parlay::sequence<uint> SA = suffixArrayRec(ss, k, findLCPs);
  //cout << "here2" << endl;
  return SA;
}

