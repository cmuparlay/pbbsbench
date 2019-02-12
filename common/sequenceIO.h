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

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "IO.h"
#include "../pbbslib/sequence_ops.h"

namespace benchIO {
  using namespace std;
  using namespace pbbs;
  
  typedef pair<int,int> intPair;
  typedef pair<unsigned int, unsigned int> uintPair;
  typedef pair<unsigned int, int> uintIntPair;
  typedef pair<long,long> longPair;
  typedef pair<char*,long> stringIntPair;

  enum elementType { none, intType, intPairT, stringIntPairT, doubleT, stringT};
  //elementType dataType(long a) { return longT;}
  elementType dataType(long a) { return intType;}
  elementType dataType(int a) { return intType;}
  elementType dataType(unsigned int a) { return intType;}
  elementType dataType(double a) { return doubleT;}
  elementType dataType(char* a) { return stringT;}
  elementType dataType(intPair a) { return intPairT;}
  elementType dataType(uintPair a) { return intPairT;}
  elementType dataType(uintIntPair a) { return intPairT;}
  elementType dataType(longPair a) { return intPairT;}
  elementType dataType(stringIntPair a) { return stringIntPairT;}

  string seqHeader(elementType dt) {
    switch (dt) {
    case intType: return "sequenceInt";
    case doubleT: return "sequenceDouble";
    case stringT: return "sequenceChar";
    case intPairT: return "sequenceIntPair";
    case stringIntPairT: return "sequenceStringIntPair";
    default: 
      cout << "writeArrayToFile: type not supported" << endl; 
      abort();
    }
  }

  elementType elementTypeFromString(string s) {
    if (s == "double") return doubleT;
    else if (s == "string") return stringT;
    else if (s == "int") return intType;
    else return none;
  }

  struct seqData { 
    void* A;
    long n;
    elementType dt; 
    char* O; // used for strings to store pointer to character array
    seqData(void* _A, long _n, elementType _dt) : 
      A(_A), O(NULL), n(_n), dt(_dt) {}
    seqData(void* _A, char* _O, long _n, elementType _dt) : 
      A(_A), O(_O), n(_n), dt(_dt) {}
    void del() {
      if (O) free(O);
      free(A);
    }
  };

  seqData readSequenceFromFile(char* fileName) {
    sequence<char> S = readStringFromFile(fileName);
    sequence<char*> W = stringToWords(S);
    char* header = W[0];
    long n = W.size()-1;

    if (header == seqHeader(intType)) {
      int* A = pbbs::new_array_no_init<int>(n);
      parallel_for(0, n, [&] (long i) {
	  A[i] = atoi(W[i+1]);});
      return seqData((void*) A, n, intType);
    } else if (header == seqHeader(doubleT)) {
      double* A = pbbs::new_array_no_init<double>(n);
      parallel_for(0, n, [&] (long i) {
	  A[i] = atof(W[i+1]);});
      return seqData((void*) A, n, doubleT);
    } else if (header == seqHeader(stringT)) {
      char** A = pbbs::new_array_no_init<char*>(n);
      parallel_for(0, n, [&] (long i) {
	  A[i] = W[i+1];});
      char* s = S.to_array();
      return seqData((void*) A, s, n, stringT);
    } else if (header == seqHeader(intPairT)) {
      n = n/2;
      intPair* A = pbbs::new_array_no_init<intPair>(n);
      parallel_for (0, n, [&] (long i) {
	A[i].first = atoi(W[2*i+1]);
	A[i].second = atoi(W[2*i+2]);
	});
      return seqData((void*) A, n, intPairT);
    } else if (header == seqHeader(stringIntPairT)) {
      n = n/2;
      stringIntPair* A = pbbs::new_array_no_init<stringIntPair>(n);
      parallel_for (0, n, [&] (long i) {
	A[i].first = W[2*i+1];
	A[i].second = atoi(W[2*i+2]);
	});
      char* s = S.to_array();
      return seqData((void*) A, s, n, stringIntPairT);
    }
    abort();
  }

  template <class T>
  int writeSequenceToFile(sequence<T> A, char* fileName) {
    elementType tp = dataType(A[0]);
    return writeArrayToFile(seqHeader(tp), A.begin(), A.size(), fileName);
  }

};
