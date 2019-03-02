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
#include "../pbbslib/sequence_ops.h"
#include "../pbbslib/parallel.h"
#include "../pbbslib/get_time.h"

namespace benchIO {
  using namespace std;
  using namespace pbbs;

  inline bool isSpace(char c) {
    switch (c)  {
    case '\r': 
    case '\t': 
    case '\n': 
    case 0:
    case ' ' : return true;
    default : return false;
    }
  }

  // parallel code for converting a string to word pointers
  // side effects string by setting to null after each word
  template <class Seq>
  sequence<char*> stringToWords(Seq &Str) {
    size_t n = Str.size();
    
    parallel_for(0, n, [&] (long i) {
	if (isSpace(Str[i])) Str[i] = 0;}); 

    // mark start of words
    sequence<bool> FL(n, [&] (long i) {
	return (i==0) ? Str[0] : Str[i] && !Str[i-1];});
    
    // offset for each start of word
    sequence<long> Offsets = pbbs::pack_index<long>(FL);

    // pointer to each start of word
    sequence<char*> SA(Offsets.size(), [&] (long j) {
	return Str.begin() + Offsets[j];});
    
    return SA;
  }

  int writeStringToFile(char* S, long n, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file.write(S, n);
    file.close();
    return 0;
  }

  inline int xToStringLen(long a) { return 21;}
  inline void xToString(char* s, long a) { sprintf(s,"%ld",a);}

  inline int xToStringLen(unsigned long a) { return 21;}
  inline void xToString(char* s, unsigned long a) { sprintf(s,"%lu",a);}

  inline uint xToStringLen(uint a) { return 12;}
  inline void xToString(char* s, uint a) { sprintf(s,"%u",a);}

  inline int xToStringLen(int a) { return 12;}
  inline void xToString(char* s, int a) { sprintf(s,"%d",a);}

  inline int xToStringLen(double a) { return 18;}
  inline void xToString(char* s, double a) { sprintf(s,"%.11le", a);}

  inline int xToStringLen(char* a) { return strlen(a)+1;}
  inline void xToString(char* s, char* a) { sprintf(s,"%s",a);}

  template <class A, class B>
  inline int xToStringLen(pair<A,B> a) { 
    return xToStringLen(a.first) + xToStringLen(a.second) + 1;
  }

  template <class A, class B>
  inline void xToString(char* s, pair<A,B> a) { 
    int l = xToStringLen(a.first);
    xToString(s, a.first);
    s[l] = ' ';
    xToString(s+l+1, a.second);
  }

  template <class Seq>
  sequence<char> seqToString(Seq const &A) {
    size_t n = A.size();
    sequence<long> L(n, [&] (size_t i) -> size_t {
	return xToStringLen(A[i])+1;});
    size_t m;
    std::tie(L,m) = pbbs::scan(std::move(L), addm<size_t>());

    sequence<char> B(m+1, (char) 0);
    char* Bs = B.begin();

    parallel_for(0, n-1, [&] (long i) {
      xToString(Bs + L[i], A[i]);
      Bs[L[i+1] - 1] = '\n';
      });
    xToString(Bs + L[n-1], A[n-1]);
    Bs[m] = Bs[m-1] = '\n';
    
    sequence<char> C = pbbs::filter(B, [&] (char c) {return c != 0;}); 
    C[C.size()-1] = 0;
    return C;
  }

  template <class T>
  void writeSeqToStream(ofstream& os, sequence<T> const &A) {
    size_t bsize = 10000000;
    size_t offset = 0;
    size_t n = A.size();
    while (offset < n) {
      // Generates a string for a sequence of size at most bsize
      // and then wrties it to the output stream
      sequence<char> S = seqToString(A.slice(offset, min(offset + bsize, n)));
      os.write(S.begin(), S.size()-1);
      offset += bsize;
    }
  }

  template <class T>
  int writeSeqToFile(string header, sequence<T> const &A, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file << header << endl;
    writeSeqToStream(file, A);
    file.close();
    return 0;
  }

  sequence<char> readStringFromFile(char *fileName) {
    ifstream file (fileName, ios::in | ios::binary | ios::ate);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      abort();
    }
    long end = file.tellg();
    file.seekg (0, ios::beg);
    long n = end - file.tellg();
    sequence<char> bytes(n, (char) 0);
    file.read (bytes.begin(), n);
    file.close();
    return bytes;
  }

  string intHeaderIO = "sequenceInt";

  template <class T>
  int writeIntSeqToFile(sequence<T> const &A, char* fileName) {
    return writeSeqToFile(intHeaderIO, A, fileName);
  }

  template <class T>
  sequence<T> readIntSeqFromFile(char *fileName) {
    sequence<char> S = readStringFromFile(fileName);
    sequence<char*> W = stringToWords(S);
    string header = (string) W[0];
    if (header != intHeaderIO) {
      cout << "readIntSeqFromFile: bad input" << endl;
      abort();
    }
    long n = W.size()-1;
    sequence<T> A(n, [&] (long i) {return atol(W[i+1]);});
    return A;
  }
};

