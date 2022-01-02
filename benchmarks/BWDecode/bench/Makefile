include common/parallelDefs

BNCHMRK = bw

CHECKFILES = $(BNCHMRK)Check.o

COMMON = 

INCLUDE = 

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(BNCHMRK)Check : $(CHECKFILES)
	$(CC) $(LFLAGS) -o $@ $(CHECKFILES)

clean :
	rm -f $(BNCHMRK)Check *.o
