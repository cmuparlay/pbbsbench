
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

  edgeArray<intT> G = readEdgeArrayFromFile<intT>(iFile);
  edge<intT>* E = G.E;
  parallel_for(long i=0;i<G.nonZeros;i++) {
    if(E[i].u > E[i].v) swap(E[i].u,E[i].v);
  }
  edge<intT> *EE = newA(edge<intT>,G.nonZeros);
  intT mm = sequence::filter(E,EE,G.nonZeros,nEQF<intT>());
  free(E);
  G.E = EE; G.nonZeros = mm;
  edgeArray<intT> F = remDuplicates2(G);
  G.del();

  writeEdgeArrayToFile<intT>(F,oFile);
}
