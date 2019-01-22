#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "parallel.h"
using namespace benchIO;
using namespace std;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile> inFile permutation");
  char* iFile = P.getArgument(1);
  char* iFile2 = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<uintT> G = readGraphFromFile<uintT>(iFile);
  _seq<uintT> I = readIntArrayFromFile<uintT>(iFile2);
  graph<uintT> GG = graphReorder(G,I.A);
  I.del();
  writeGraphToFile<uintT>(GG,oFile);
  GG.del();
}
