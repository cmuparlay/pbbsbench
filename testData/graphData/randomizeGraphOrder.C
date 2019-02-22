
#include "pbbslib/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile> inFile permutation");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<uintT> G = readGraphFromFile<uintT>(iFile);
  graph<uintT> GG = graphReorder(G,(uintT*)NULL);
  writeGraphToFile<uintT>(GG,oFile);
  GG.del();
}
