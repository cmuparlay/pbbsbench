#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
using namespace benchIO;

sequence<char*> trigramWords(size_t s, size_t e);

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<size> <outfile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  sequence<char*> A = trigramWords(0,in.first);
  return writeSequenceToFile(A, in.second);
}
