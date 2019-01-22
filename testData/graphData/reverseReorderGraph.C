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
  uintT* reverseI = newA(uintT,I.n);
  parallel_for(long i=0;i<I.n;i++) reverseI[I.A[i]] = i;
  I.del();
  graph<uintT> GG = graphReorder(G,reverseI);
  free(reverseI);
  writeGraphToFile<uintT>(GG,oFile);
  GG.del();
}
