// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2017 Guy Blelloch and the PBBS team
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
#include <algorithm>
#include <math.h>
#include "get_time.h"
#include "sequence_ops.h"
#include "integer_sort.h"
#include "quicksort.h"
#include "sample_sort.h"
#include "parallel.h"
#include "SA.h"
using namespace std;

#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))

//#define printInfo

#ifdef printInfo
#define nextTimeM(_str) nextTime(_str)
#else
#define nextTimeM(_str) 
#endif

typedef int intT;
typedef unsigned int uintT;
typedef unsigned char uchar;
typedef pair<uintT,uintT> intpair;
typedef unsigned __int128 long_int;

auto maxuint = [] (uintT a, uintT b) {return std::max<uintT>(a,b);};
    
struct seg {
  uintT start;
  uintT length;
  seg(uintT s, uintT l) : start(s), length(l) {}
};

struct isSeg {bool operator() (seg s) {return s.length > 1;}};

struct pairCompF {
  bool operator() (intpair A, intpair B) { return A.first < B.first;}};

template <typename kvpair>
void splitSegment(seg *segOut, uintT start, uintT l, uintT* ranks, kvpair *Cs) {
  //void splitSegment(seg *segOut, uintT start, uintT l, uintT* ranks, intpair *Cs) {
  if (l < 5000) { // sequential version

    uintT name = 0;
    ranks[Cs[0].second] = name + start + 1;
    for (uintT i=1; i < l; i++) {
      if (Cs[i-1].first != Cs[i].first) name = i;
      ranks[Cs[i].second] = name + start + 1;
    }

    name = 0;
    for (uintT i=1; i < l; i++) {
      if (Cs[i-1].first != Cs[i].first) {
	segOut[i-1] = seg(name+start,i-name);
	name = i;
      } else segOut[i-1] = seg(0,0);
    }
    segOut[l-1] = seg(name+start,l-name);

  } else { // parallel version
    uintT *names = newA(uintT,l);

    // mark start of each segment with equal keys
    parallel_for (1, l, [&] (size_t i) {
	names[i] = (Cs[i].first != Cs[i-1].first) ? i : 0;});
    names[0] = 0;

    // scan start i across each segment
    sequence<uintT> snames(names,l);
    pbbs::scan(snames, snames, maxuint, (uintT) 0,
	       pbbs::fl_scan_inclusive);

    // write new rank into original location
    parallel_for (0, l, [&] (size_t i) {
	ranks[Cs[i].second] = names[i]+start+1;});

    // get starts and lengths of new segments
    parallel_for (1, l, [&] (size_t i) {
	if (names[i] == i) 
	  segOut[i-1] = seg(start+names[i-1],i-names[i-1]);
	else segOut[i-1] = seg(0,0);
      });
    segOut[l-1] = seg(start+names[l-1],l-names[l-1]);

    free(names);
  }
}  

intpair* splitSegmentTop(seg *segOut, uintT n, uintT* ranks, long_int *Cs) {
  uintT *names = newA(uintT,n);
  size_t mask = ((((size_t) 1) << 32) - 1);

  // mark start of each segment with equal keys
  parallel_for (1, n, [&] (size_t i) {
      names[i] = ((Cs[i] >> 32) != (Cs[i-1] >> 32)) ? i : 0;});
  names[0] = 0;
  nextTimeM("names");

  // scan start i across each segment
  sequence<uintT> snames(names,n);
  pbbs::scan(snames, snames, maxuint, (uintT) 0,
	     pbbs::fl_scan_inclusive);

  intpair *C = newA(intpair,n);
  // write new rank into original location
  parallel_for (0, n, [&] (size_t i) {
      ranks[Cs[i] & mask] = names[i]+1;
      C[i].second = Cs[i] & mask;
    });
  nextTimeM("write rank and copy");

  // get starts and lengths of new segments
  parallel_for (1, n, [&] (size_t i) { 
      if (names[i] == i) 
	segOut[i-1] = seg(names[i-1],i-names[i-1]);
      else segOut[i-1] = seg(0,0);
    });
  segOut[n-1] = seg(names[n-1],n-names[n-1]);

  nextTimeM("segments");
  free(names);
  return C;
}

uintT* suffixArrayInternal(uchar* ss, long n) { 
  startTime();

  // renumber characters densely
  // start numbering at 1 leaving 0 to indicate end-of-string
  // pad the end of string with 0s
  size_t pad = 48;
  uchar *s = newA(uchar,n+pad);
  uintT flags[256];
  for (uintT i=0; i < 256; i++) flags[i] = 0;
  parallel_for (0, n, [&] (size_t i) {
      if (!flags[ss[i]]) flags[ss[i]] = 1;});
  sequence<uintT> sflags(flags,256);
  auto add = [&] (uintT a, uintT b) {return a + b;};
  uintT m = pbbs::scan(sflags, sflags, add, (uint) 1);
#ifdef printInfo
  cout << "m = " << m << " n = " << n << endl;
#endif
  parallel_for (0, n, [&] (size_t i) {
      s[i] = flags[ss[i]];});
  for (uintT i=n; i < n + pad; i++) s[i] = 0;

  // pack characters into 128-bit word, along with the location i
  // 96 bits for characters, and 32 for location
  double logm = log2((double) m);
  uintT nchars = floor(96.0/logm); 

  //long_int *Cl = newA(long_int,n);

  sequence<long_int> Cl(n, [&] (size_t i) {
      long_int r = s[i];
      for (uintT j=1; j < nchars; j++) r = r*m + s[i+j];
      return (r << 32) + i;
    });
  free(s);
  nextTimeM("copy");

  // sort based on packed words
  pbbs::sample_sort(Cl, std::less<long_int>(), true);
  nextTimeM("sort");

  // identify segments of equal values
  uintT *ranks = newA(uintT,n);
  seg *segOuts = newA(seg,n);
  intpair *C = splitSegmentTop(segOuts, n, ranks, Cl.start());
  Cl.clear();
  nextTimeM("split");

  uintT *offsets = newA(uintT,n);
  seg *segments= newA(seg,n);

  uintT offset = nchars;
  uint round =0;
  uintT nKeys = n;
  while (1) {
    if (round++ > 40) {
      cout << "Suffix Array:  Too many rounds" << endl;
      abort();
    }
    
    sequence<seg> Segs = pbbs::filter(sequence<seg>(segOuts,nKeys),isSeg(),
				      pbbs::no_flag, segments);
    uintT nSegs = Segs.size();
    if (nSegs == 0) break;
    
#ifdef printInfo
    cout << "nSegs = " << nSegs << " nKeys = " << nKeys 
	 << " common length = " << offset << endl;
#endif
    nextTimeM("filter and scan");    

    parallel_for (0, nSegs, [&] (size_t i) {
	uintT start = segments[i].start;
	intpair *Ci = C + start;
	uintT l = segments[i].length;
	offsets[i] = l;
	parallel_for (0, l, [&] (size_t j) {
	    uintT o = Ci[j].second+offset;
	    Ci[j].first = (o >= n) ? 0 : ranks[o]; 
	  }, 100);
	auto less = [&] (intpair A, intpair B) {return A.first < B.first;};
	if (l >= n/10) pbbs::sample_sort(Ci, l, less);
	else pbbs::quicksort(Ci, l, less);
      }, 100);
    nextTimeM("sort");

    sequence<uintT> soffsets(offsets, nSegs);
    nKeys = pbbs::scan_add(soffsets, soffsets);

    parallel_for (0, nSegs, [&] (size_t i) {
	uintT start = segments[i].start;
	splitSegment(segOuts + offsets[i], start, segments[i].length, 
		     ranks, C + start);
      }, 100);
    nextTimeM("split");

    offset = 2 * offset;
  }
  parallel_for (0, n, [&] (size_t i) {ranks[i] = C[i].second;});
  free(C); free(segOuts); free(segments); free(offsets); 
  return ranks;
}

sequence<uintT> suffixArray(sequence<unsigned char> s) {
  //sequence<uchar> ss(s.size()+1, [&] (size_t i) {return (i==s.size()) ? 0 : s[i];});
  return sequence<uintT>(suffixArrayInternal(s.start(), s.size()), s.size());
}
