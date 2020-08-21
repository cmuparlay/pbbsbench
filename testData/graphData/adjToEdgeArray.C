#include "common/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile> <inflile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<size_t> G = readGraphFromFile<size_t>(iFile);
  edgeArray<size_t> EA = edgesFromGraph(G);
  EA = makeSymmetric(EA);
  writeEdgeArrayToFile(EA,oFile);
}
