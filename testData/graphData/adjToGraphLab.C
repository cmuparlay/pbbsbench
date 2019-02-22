
#include "pbbslib/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
#include "pbbslib/parallel.h"
#include <iostream>
#include <sstream>
using namespace benchIO;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<intT> G = readGraphFromFile<intT>(iFile);
  edgeArray<intT> E = edgesFromGraph<intT>(G);
  G.del();
  G = graphFromEdges<intT>(E,1);
  E.del();

  ofstream out(oFile);

  uintT n = G.n;
  for(uintT i=0;i<n;i++) {
    stringstream ss; 
    uintT d = 0;
    for(uintT j=0;j<G.V[i].degree;j++)
      if (G.V[i].Neighbors[j] > i) d++;
    ss << i << " " << d << " ";
    for(uintT j = 0;j<G.V[i].degree;j++){
      if(G.V[i].Neighbors[j] > i) ss << G.V[i].Neighbors[j] << " ";
    }
    ss << "\n";
    out << ss.str();
  }

  G.del();
  out.close();
}
