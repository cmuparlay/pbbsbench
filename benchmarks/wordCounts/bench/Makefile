include common/parallelDefs

BNCHMRK = wc

CHECKFILES = serialwc.o $(BNCHMRK)Check.o

COMMON = 

INCLUDE = 

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(BNCHMRK)Check : $(CHECKFILES)
	$(CC) -o $@ $(CHECKFILES) $(LFLAGS)

clean :
	rm -f wcCheck *.o
