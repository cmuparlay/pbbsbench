include common/seqDefs
LFLAGS := $(LFLAGS) -fopenmp
CFLAGS := $(CFLAGS) -fopenmp -DPARLAY_SEQUENTIAL -DNDEBUG
BENCH = sort

include common/MakeBench
