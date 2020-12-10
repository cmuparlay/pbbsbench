#include "sequenceData.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
#include "parlay/random.h"
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
  parlay::random rnd(0);
  parlay::sequence<long> A = parlay::tabulate(n, [&] (size_t i) -> long {
      size_t x = rnd.ith_rand(i);
      return (x % 2 == 0) ? val : x % r;});
  
  switch(tp) {
  case intType: {
    return writeSequenceToFile(A, fname);
  }
  case doubleT: {
    return writeSequenceToFile(parlay::map(A, [] (size_t i) {return (double) i;}), fname);
  }
  default: cout << "genSeqRand: not a valid type" << endl;
  }
}
