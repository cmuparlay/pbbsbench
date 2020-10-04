#include "parlay/parallel.h"
#include "common/parse_command_line.h"
#include "parlay/io.h"
#include "common/sequenceIO.h"
using namespace benchIO;

char* trigramString(size_t s, size_t e);

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<size> <outfile>");
  pair<size_t, char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* S = trigramString(0, n);
  parlay::chars_to_file(parlay::sequence<char>(S,S+n), in.second);
  return 0;
}
