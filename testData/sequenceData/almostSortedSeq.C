#include "sequenceData.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
using namespace dataGen;
using namespace benchIO;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-r <swaps>] [-t {int,double}] <size> <outfile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  size_t n = in.first;
  char* fname = in.second;
  size_t swaps = P.getOptionIntValue("-r",floor(sqrt((float) n)));
  
  switch(tp) {
  case intType: 
    return writeSequenceToFile(almostSorted<long>(0, n, swaps), fname);
  case doubleT: 
    return writeSequenceToFile(almostSorted<double>(0, n, swaps), fname);
  default:
    cout << "genSeqRand: not a valid type" << endl;
  }
}
