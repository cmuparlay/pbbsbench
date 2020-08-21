#include "sequenceData.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
using namespace dataGen;
using namespace benchIO;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-t {int,double}] <size> <outfile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  switch(tp) {
  case intType: 
    return writeSequenceToFile(expDist<long>(0,n), fname);
  case doubleT: 
    return writeSequenceToFile(expDist<double>(0,n), fname);
  default:
    cout << "genSeqRand: not a valid type" << endl;
  }
}
