
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
  commandLine P(argc,argv,"<inFile> <idxFile> <adjFile> <configFile>");
  char* iFile = P.getArgument(3);
  char* idxFile = P.getArgument(2);
  char* adjFile = P.getArgument(1);
  char* configFile = P.getArgument(0);

  wghGraph<uintT> G = readWghGraphFromFile<uintT>(iFile);

  ofstream idx(idxFile, ofstream::out | ios::binary);
  ofstream adj(adjFile, ofstream::out | ios::binary);
  ofstream config(configFile, ofstream::out);
  config << G.n;
  config.close();
  uintT* In = G.allocatedInplace;
  uintT* offsets = In+2;
  uintT* edges = In+2+G.n;
  
  idx.write((char*)offsets,sizeof(uintT)*G.n);
  adj.write((char*)edges,sizeof(uintT)*2*G.m); //edges and weights
  idx.close();
  adj.close();
}
