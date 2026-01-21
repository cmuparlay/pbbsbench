include common/parallelDefs

INCLUDE = -Icommon

SACheck: SACheck.C
	$(CC) $(CFLAGS) $(INCLUDE) -o SACheck SACheck.C $(LFLAGS)

SATime.o: SATime.C
	$(CC) $(CFLAGS) $(INCLUDE) -o SATime.o -c SATime.C

clean :
	rm -f SACheck *.o
