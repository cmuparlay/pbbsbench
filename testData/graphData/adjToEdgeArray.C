#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "parallel.h"
using namespace benchIO;
using namespace std;

int parallel_main(int argc, char* argv[]) {
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
