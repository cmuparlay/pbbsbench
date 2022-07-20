ifeq (, $(shell which jemalloc-config))
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -L$(JEMALLOCLD) -ljemalloc 
endif

CC = g++
CFLAGS = -DPARLAY_TBB -mcx16 -O3 -std=c++17 -DPARLAY_USE_STD_ALLOC
LFLAGS = $(JEMALLOC)

BENCH = neighbors
TIME = ../bench/$(BENCH)Time.C
CHECK = $(BENCH)Check
INCLUDE = -Icommon
CGALINCLUDE = -I/home/magdalen/CGAL-5.2/include -I/home/magdalen/tbb/include -DCGAL_LINKED_WITH_TBB
CGALLIB = /home/magdalen/tbb/build/libtbb.so

all : $(BENCH) testInputs
	cd ../bench; make -s $(CHECK)

$(BENCH) : $(TIME) $(BENCH).h $(REQUIRE)
	$(CC) $(CFLAGS) $(INCLUDE) $(CGALINCLUDE) -include $(BENCH).h -o $(BENCH).o -c $(TIME)
	$(CC) -o $(BENCH) $(BENCH).o $(CGALLIB) $(LFLAGS)

test : test.cpp
	g++ -O3 $(CGALINCLUDE) -o test.o -c test.cpp
	g++ -o test test.o $(CGALLIB)

testInputs : ../bench/testInputs
	cp ../bench/testInputs .

clean :
	rm -f $(BENCH) *.o test

cleanall : clean
	rm -f testInputs; cd ../bench; make -s clean

