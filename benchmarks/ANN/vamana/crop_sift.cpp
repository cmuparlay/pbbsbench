#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "../utils/types.h"
#include "../utils/NSGDist.h"
#include "../utils/parse_files.h"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "usage: crop_sift <base> <oF>" << std::endl;
    return 1;
  }
  auto [fileptr, length] = mmapStringFromFile(argv[1]);
  parlay::sequence<int> preamble = {1000000, 128};

  int* pr = preamble.begin();
  uint8_t* data = (uint8_t*)(fileptr+8);
  std::ofstream writer;
  writer.open(argv[2], std::ios::binary | std::ios::out);
  writer.write((char *) pr, 2*sizeof(int));
  writer.write((char *) data, 1000000*128*sizeof(uint8_t));
  writer.close();

  return 0;
}