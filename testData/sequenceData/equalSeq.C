#include "sequenceData.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
using namespace benchIO;
using namespace dataGen;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-r <range>] [-t {int,double}] <size> <outfile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  size_t n = in.first;
  char* fname = in.second;
  size_t val = 1723451871;
  size_t r = P.getOptionLongValue("-r", 0);
  if (r > 0) val = val % r;
  sequence<long> A;
  sequence<double> B;
  
  switch(tp) {
  case intType: {
    auto A = parlay::tabulate(n, [&] (size_t i) -> long {return val;});
    return writeSequenceToFile(A, fname);
  }
  case doubleT: {
    auto B = parlay::tabulate(n, [&] (size_t i) -> double {return (double) val;});
    return writeSequenceToFile(B, fname);
  }
  default: cout << "genSeqRand: not a valid type" << endl;
  }
}
