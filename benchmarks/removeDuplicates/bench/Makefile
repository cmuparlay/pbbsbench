include common/parallelDefs

BNCHMRK = dedup

$(BNCHMRK)Check : $(BNCHMRK)Check.C
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $(BNCHMRK)Check.C

clean :
	rm -f $(BNCHMRK)Check *.o
