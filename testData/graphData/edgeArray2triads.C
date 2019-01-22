#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "graphUtils.h"
#include "parallel.h"
using namespace benchIO;
using namespace std;

template <class intT>
struct directed {bool operator() (edge<intT> e) {return (e.u < e.v);}};


int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-o <outFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  edgeArray<intT> EA = readEdgeArrayFromFile<intT>(iFile);
  cout << "read edge array\n";
  edgeArray<intT> EA2 = makeSymmetric(EA);
  EA.del();
  cout << "made symmetric\n";
  edge<intT>* E = newA(edge<intT>,EA2.nonZeros);
  intT new_m = sequence::filter(EA2.E,E,EA2.nonZeros,directed<intT>());
  cout << "filtered\n";
  stringstream ss;
  ss << max(EA2.numRows,EA2.numCols);
  cout << "writing to file...\n";
  EA2.del();

  writeArrayToFile(ss.str(), E, new_m, oFile);
  free(E);
}
