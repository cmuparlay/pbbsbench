include common/parallelDefsANN

REQUIRE = ../utils/beamSearch.h 
BENCH = range

include common/MakeBench

compute_range_groundtruth : compute_range_groundtruth.cpp
	$(CC) $(CFLAGS) -o compute_range_groundtruth compute_range_groundtruth.cpp $(LFLAGS) 