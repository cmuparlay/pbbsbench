include common/parallelDefs

INCLUDE = -Icommon

SACheck: SACheck.C 
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDE) -o SACheck SACheck.C

SATime.o: SATime.C
	$(CC) $(CFLAGS) $(INCLUDE) -o SATime.o -c SATime.C

clean :
	rm -f SACheck *.o
