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
  size_t r = P.getOptionLongValue("-r",n);

  switch(tp) {
  case intType: return writeSequenceToFile(randIntRange<long>((size_t) 0,n,r),
					   fname);
  case doubleT: return writeSequenceToFile(rand<double>((size_t) 0, n), fname);
  default: cout << "genSeqRand: not a valid type" << endl;
  }
}
