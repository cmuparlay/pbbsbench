
#include "pbbslib/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<intT> G = readGraphFromFile<intT>(iFile);

  edgeArray<intT> EA = edgesFromGraph(G);
  G.del();

  EA = makeSymmetric(EA);
  writeEdgeArrayToFile<intT>(EA,oFile);
  EA.del();
}
