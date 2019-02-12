#include "sequenceData.h"
#include "common/sequenceIO.h"
#include "pbbslib/parse_command_line.h"
using namespace benchIO;
using namespace dataGen;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-r <range>] [-t {int,double}] <size> <outfile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  size_t n = in.first;
  char* fname = in.second;
  size_t defaultVal = 1723451871;
  size_t r = P.getOptionLongValue("-v", defaultVal);
  sequence<long> A;
  sequence<double> B;
  
  switch(tp) {
  case intType:
    A = sequence<long>(n, [&] (size_t i) {return r;});
    return writeSequenceToFile(A, fname);
  case doubleT:
    B = sequence<double>(n, [&] (size_t i) {return (double) r;});
    return writeSequenceToFile(B, fname);
  default: cout << "genSeqRand: not a valid type" << endl;
  }
}
