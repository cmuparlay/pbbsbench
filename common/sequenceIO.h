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
#include "../parlay/primitives.h"
#include "../parlay/parallel_io.h"

namespace benchIO {
  using namespace std;
  using parlay::sequence;
  using parlay::tabulate;
  using parlay::make_slice;

  typedef unsigned int uint;
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
  elementType dataType(uint a) { return intType;}
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

  template <typename Range>
  elementType elementTypeFromHeader(Range R) {
    string s(R.begin(), R.end());
    if (s == "sequenceInt") return intType;
    else if (s == "sequenceDouble") return doubleT;
    else if (s == "sequenceChar") return stringT;
    else if (s == "sequenceIntPair") return intPairT;
    else if (s == "sequenceStringIntPair") return stringIntPairT;
    else if (s == "sequenceDoublePair") return doublePairT;
    else return none;
  }

  template <typename Range>
  elementType elementTypeFromString(Range R) {
    string s(R.begin(), R.end());
    if (s == "double") return doubleT;
    else if (s == "string") return stringT;
    else if (s == "int") return intType;
    else return none;
  }

  long read_long(sequence<char> S) {
    return parlay::char_range_to_l(S);}

  double read_double(sequence<char> S) {
    return parlay::char_range_to_d(S);}

  using charseq_slice = decltype(make_slice(sequence<sequence<char>>()));
  
  template <typename T>
  sequence<T> parseElements(charseq_slice S);

  // specialized parsing functions
  template<>
  sequence<double> parseElements<double>(charseq_slice S) {
    return tabulate(S.size(), [&] (long i) -> double {return read_double(S[i]);});
  }

  template<>
  sequence<int> parseElements<int>(charseq_slice S) {
    return tabulate(S.size(), [&] (long i) -> int {return (int) read_long(S[i]);});
  }

  template<>
  sequence<uint> parseElements<uint>(charseq_slice S) {
    return tabulate(S.size(), [&] (long i) -> uint {return (uint) read_long(S[i]);});
  }

  template<>
  sequence<intPair> parseElements<intPair>(charseq_slice S) {
    return tabulate((S.size())/2, [&] (long i) -> intPair {
      return std::make_pair((int) read_long(S[2*i]), (int) read_long(S[2*i+1]));});
  }

  template<>
  sequence<uintPair> parseElements<uintPair>(charseq_slice S) {
    return tabulate((S.size())/2, [&] (long i) -> uintPair {
      return std::make_pair((uint) read_long(S[2*i]), (uint) read_long(S[2*i+1]));});
  }

  template<>
  sequence<doublePair> parseElements<doublePair>(charseq_slice S) {
    return tabulate((S.size())/2, [&] (long i) -> doublePair {
      return std::make_pair(read_double(S[2*i]), read_double(S[2*i+1]));});
  }

  template<>
  sequence<char*> parseElements<char*>(charseq_slice S) {
    return sequence<char*>(0);
  }

  template<>
  sequence<stringIntPair> parseElements<stringIntPair>(charseq_slice S) {
    return sequence<stringIntPair>(0);
  }  

  sequence<sequence<char>> get_tokens(char const *fileName) {
    sequence<char> S = parlay::char_seq_from_file(fileName);
    return parlay::tokens(S, benchIO::is_space);
  }

  template <typename T, typename CharRange>
  void check_header(CharRange& S) {
    T a;
    string header(S[0].begin(), S[0].end());
    string type_str = seqHeader(dataType(a));
    if (header != type_str) {
      cout << "bad header: expected " << type_str << " got " << header << endl;
      abort();
    }
  }

  // reads file, tokenizes and then dispatches to specialized parsing function
  template <typename T>
  sequence<T> readSequenceFromFile(char const *fileName) {
    auto S = get_tokens(fileName);
    check_header<T>(S[0]);
    return parseElements<T>(S.cut(1,S.size()));
  }
  
  template <class T>
  int writeSequenceToFile(sequence<T> const &A, char const *fileName) {
    elementType tp = dataType(A[0]);
    return writeSeqToFile(seqHeader(tp), A, fileName);
  }

};
