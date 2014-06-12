include ../../Makefile.inc

GCC      = gcc
GCCFLAGS = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wa,-q

all: ../../lib/libtpl.a

../../lib/libtpl.a: tpl.o
	ar cr ../../lib/libtpl.a tpl.o

tpl.o: tpl.c tpl.h
	$(GCC) $(GCCFLAGS) tpl.c -o tpl.o

clean:
	rm -f tpl.o
	rm -f ../../lib/libtpl.a
