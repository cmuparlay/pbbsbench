#include "pbbslib/parallel.h"
#include "pbbslib/parse_command_line.h"
#include "common/sequenceIO.h"
using namespace benchIO;

char* trigramString(size_t s, size_t e);

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<size> <outfile>");
  pair<size_t, char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* S = trigramString(0, n);
  return writeStringToFile(S, n, in.second);
}
