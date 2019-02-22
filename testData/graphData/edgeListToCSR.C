#include <iostream>
#include <algorithm>
#include "utils.h"
#include "common/graph.h"
#include "pbbslib/parallel.h"

#include "common/graphIO.h"
#include "pbbslib/parse_command_line.h"
#include "common/graphUtils.h"
using namespace std;
using namespace benchIO;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  edgeArray<intT> EA = readEdgeArrayFromFile<intT>(iFile);
  graph<intT> G = graphFromEdges(EA,1);
  writeGraphToFile(G, oFile);
}
