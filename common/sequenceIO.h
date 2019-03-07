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
  typedef pair<double,double> doublePair;

  enum elementType { none, intType, intPairT, doublePairT,
		     stringIntPairT, doubleT, stringT};
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
  elementType dataType(doublePair a) { return doublePairT;}

  string seqHeader(elementType dt) {
    switch (dt) {
    case intType: return "sequenceInt";
    case doubleT: return "sequenceDouble";
    case stringT: return "sequenceChar";
    case intPairT: return "sequenceIntPair";
    case stringIntPairT: return "sequenceStringIntPair";
    case doublePairT: return "sequenceDoublePair";
    default: 
      cout << "writeSeqToFile: type not supported" << endl; 
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
    sequence<char> s;  // backup store if A has pointers to strings
    seqData(void* A, long n, elementType dt) 
      : A(A), s(sequence<char>()), n(n), dt(dt) {}
    seqData(void* A, sequence<char> s, long n, elementType dt)
      : A(A), s(std::move(s)), n(n), dt(dt) {}
    void clear() {pbbs:free_array(A); A = NULL; n = 0;}
  };

  seqData readSequenceFromFile(char* fileName) {
    sequence<char> S = readStringFromFile(fileName);
    sequence<char*> W = stringToWords(S);
    char* header = W[0];
    long n = W.size()-1;

    if (header == seqHeader(intType)) {
      sequence<int> A(n, [&] (long i) {return atoi(W[i+1]);});
      return seqData((void*) A.to_array(), n, intType);
    } else if (header == seqHeader(doubleT)) {
      sequence<double> A(n, [&] (long i) {return atof(W[i+1]);});
      return seqData((void*) A.to_array(), n, doubleT);
    } else if (header == seqHeader(stringT)) {
      sequence<char*> A(n, [&] (long i) {return W[i+1];});
      return seqData((void*) A.to_array(), std::move(S), n, stringT);
    } else if (header == seqHeader(intPairT)) {
      n = n/2;
      sequence<intPair> A(n, [&] (long i) {
	  return std::make_pair(atoi(W[2*i+1]), atoi(W[2*i+2]));});
      return seqData((void*) A.to_array(), n, intPairT);
    } else if (header == seqHeader(doublePairT)) {
      n = n/2;
      sequence<doublePair> A(n, [&] (long i) {
	  return std::make_pair(atof(W[2*i+1]), atof(W[2*i+2]));});
      return seqData((void*) A.to_array(), n, doublePairT);
    } else if (header == seqHeader(stringIntPairT)) {
      n = n/2;
      sequence<stringIntPair> A(n, [&] (long i) {
	  return std::make_pair(W[2*i+1], atoi(W[2*i+2]));});
      return seqData((void*) A.to_array(), std::move(S), n, stringT);
    }
    abort();
  }

  template <class T>
  int writeSequenceToFile(sequence<T> const &A, char* fileName) {
    elementType tp = dataType(A[0]);
    return writeSeqToFile(seqHeader(tp), A, fileName);
  }

};
