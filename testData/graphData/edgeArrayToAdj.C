
#include "pbbslib/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-w] -o <outFile> <inflile>");
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
