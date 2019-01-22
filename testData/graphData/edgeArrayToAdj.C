#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "parallel.h"
using namespace benchIO;
using namespace std;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-w] -o <outFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  bool weighted = P.getOptionValue("-w");
  if(!weighted) {
    edgeArray<intT> G = readEdgeArrayFromFile<intT>(iFile);
    writeGraphToFile<intT>(graphFromEdges(G,1),oFile);
  } else {
    wghEdgeArray<intT> G = readWghEdgeArrayFromFile<intT>(iFile);
    writeWghGraphToFile<intT>(wghGraphFromEdges(G),oFile);
  }
}
