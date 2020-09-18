
#include "common/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-w] -o <outFile> <inflile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  bool weighted = P.getOptionValue("-w");
  if(!weighted) {
    edgeArray<size_t> G = readEdgeArrayFromFile<size_t>(iFile);
    writeGraphToFile(graphFromEdges(G,1),oFile);
  } else {
    wghEdgeArray<size_t,double> G = readWghEdgeArrayFromFile<size_t,double>(iFile);
    writeWghGraphToFile(wghGraphFromEdges(G),oFile);
  }
}
