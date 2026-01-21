include common/parallelDefs

BNCHMRK = dedup

$(BNCHMRK)Check : $(BNCHMRK)Check.C
	$(CC) $(CFLAGS) -o $@ $(BNCHMRK)Check.C $(LFLAGS)

clean :
	rm -f $(BNCHMRK)Check *.o
