CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O0
CXXFLAGS=	-g -Wall -O0 -std=c++17
CPPFLAGS=	-DTB_DEBUG #-DGWF_DEBUG

INCLUDES=
OBJS=		kalloc.o gwf-ed.o gfa-base.o gfa-io.o gfa-sub.o
PROG=		gwf-test
LIBS=		-lz -lpthread -lm -fopenmp

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

all:$(PROG)

$(PROG):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

gfa-base.o: gfa-priv.h gfa.h kstring.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
gfa-sub.o: gfa-priv.h gfa.h kalloc.h kavl.h khash.h ksort.h
gwf-ed.o: gwfa.h kalloc.h ksort.h khashl.h kdq.h kvec.h tb.h
gwfa-lin.o: gwfa.h kalloc.h ksort.h
kalloc.o: kalloc.h
main.o: gfa.h gfa-priv.h gwfa.h ketopt.h kalloc.h kseq.h
