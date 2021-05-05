include common/seqDefs
CFLAGS = -O3 -DNDEBUG -DINPLACE -DPARLAY_OPENMP -std=c++17  -march=native  -fopenmp -pthread -latomic
LFLAGS = -ldl  /usr/lib/x86_64-linux-gnu/libtbb.so.2 -latomic -mcx16
CC = g++

BENCH = sort

include common/MakeBench
