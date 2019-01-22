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
  commandLine P(argc,argv,"<inFile> <idxFile> <adjFile> <configFile>");
  char* iFile = P.getArgument(3);
  char* idxFile = P.getArgument(2);
  char* adjFile = P.getArgument(1);
  char* configFile = P.getArgument(0);

  graph<uintT> G = readGraphFromFile<uintT>(iFile);

  ofstream idx(idxFile, ofstream::out | ios::binary);
  ofstream adj(adjFile, ofstream::out | ios::binary);
  ofstream config(configFile, ofstream::out);
  config << G.n;
  config.close();
  uintT* In = G.allocatedInplace;
  uintT* offsets = In+2;
  uintT* edges = In+2+G.n;
  idx.write((char*)offsets,sizeof(uintT)*G.n);
  adj.write((char*)edges,sizeof(uintT)*G.m);
  idx.close();
  adj.close();
}
