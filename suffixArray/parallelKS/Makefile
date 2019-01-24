include common/parallelDefs

BENCH = SA

TIME = ../bench/$(BENCH)Time
CHECK = ../bench/$(BENCH)Check
INCLUDE = -Icommon -Ipbbslib

$(BENCH) : $(BENCH).o $(TIME).o $(CHECK) testInputs
	$(CC) $(LFLAGS) $(BENCH).o $(TIME).o -o $(BENCH)

$(TIME).o : $(TIME).C
	cd ../bench; make -s SATime.o

$(CHECK) : $(CHECK).C
	cd ../bench; make -s SACheck

$(BENCH).o : $(BENCH).C
	$(CC) $(CFLAGS) $(INCLUDE) -c $(BENCH).C -o $(BENCH).o

testInputs : ../bench/testInputs
	cp ../bench/testInputs .

clean :
	rm -f $(BENCH) *.o

cleanall : clean
	rm -f testInputs

