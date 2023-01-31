include common/parallelDefsANN

REQUIRE = ../utils/beamSearch.h index.h ../utils/indexTools.h
BENCH = neighbors

include common/MakeBench

bvec_to_u8bin : bvec_to_u8bin.cpp
	$(CC) $(CFLAGS) -o bvec_to_u8bin bvec_to_u8bin.cpp $(LFLAGS) 

compute_groundtruth : compute_groundtruth.cpp
	$(CC) $(CFLAGS) -o compute_groundtruth compute_groundtruth.cpp $(LFLAGS) 

crop_sift : crop_sift.cpp
	$(CC) $(CFLAGS) -o crop_sift crop_sift.cpp $(LFLAGS) 
