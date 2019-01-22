CC = cc
CFLAGS = -O2
BINDIR = ./bin

all: ${BINDIR}/washington ${BINDIR}/ac ${BINDIR}/ak ${BINDIR}/genrmf

washington: ${BINDIR}/washington

ac: ${BINDIR}/ac

ak: ${BINDIR}/ak

genrmf: ${BINDIR}/genrmf

${BINDIR}/washington:
	cd washington; \
	${CC} washington.c -o ../${BINDIR}/washington ${CFLAGS}; \

${BINDIR}/ac:
	cd ac_dense; \
	${CC} ac.c -o ../${BINDIR}/ac ${CFLAGS}; \

${BINDIR}/ak:
	cd ak; \
	${CC} ak.c -o ../${BINDIR}/ak ${CFLAGS}; \

${BINDIR}/genrmf:
	cd genrmf; \
	make; \
	cp genrmf ../${BINDIR}

clean:
	cd ${BINDIR}; rm -f washington ac ak genrmf core; \
	cd ../genrmf; make clean
