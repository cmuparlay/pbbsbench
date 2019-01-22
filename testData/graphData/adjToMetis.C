#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "parallel.h"
#include <iostream>
#include <sstream>
using namespace benchIO;
using namespace std;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  graph<uintT> G = readGraphFromFile<uintT>(iFile);
  ofstream out(oFile);

  stringstream ss; 
  ss << G.n << " " << G.m/2 << "\n";
  out << ss.str();

  uintT n = G.n;
  for(uintT i=0;i<n;i++) {
    stringstream ss; 
    for(uintT j = 0;j<G.V[i].degree;j++){
      ss << 1+G.V[i].Neighbors[j] << " ";
    }
    ss << "\n";
    out << ss.str();
  }

  G.del();
  out.close();
}
