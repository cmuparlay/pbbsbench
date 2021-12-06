include common/parallelDefs

COMMON = common/sequenceIO.h common/IO.h common/parse_command_line.h
LIB = parlay/parallel.h
SEQUENCEGEN = $(COMMON) $(LIB) 
GENERATORS = equalSeq randomSeq almostSortedSeq almostEqualSeq exptSeq trigramSeq addDataSeq trigramString

.PHONY: all clean
all: $(GENERATORS)

classifyClean : classifyClean.C $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C 

addDataSeq : addDataSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

randomSeq : randomSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

equalSeq : equalSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

almostEqualSeq : almostEqualSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

almostSortedSeq : almostSortedSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

exptSeq : exptSeq.C sequenceData.h $(SEQUENCEGEN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $@.C

trigrams.o : trigrams.C $(SEQUENCEGEN) 
	$(CC) $(CFLAGS) -c trigrams.C

trigramSeq.o : trigramSeq.C $(SEQUENCEGEN)
	$(CC) $(CFLAGS) -c trigramSeq.C

trigramSeq : trigramSeq.o trigrams.o 
	$(CC) $(LFLAGS) -o $@ $@.o trigrams.o

trigramString.o : trigramString.C $(SEQUENCEGEN)
	$(CC) $(CFLAGS) -c trigramString.C

trigramString : trigramString.o trigrams.o 
	$(CC) $(LFLAGS) -o $@ $@.o trigrams.o

clean :
	rm -f *.o $(GENERATORS)
	make clean -s -C data
