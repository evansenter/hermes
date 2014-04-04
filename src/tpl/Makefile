all: ../../lib/libtpl.a

../../lib/libtpl.a: tpl.o
	ar cr ../../lib/libtpl.a tpl.o

tpl.o: tpl.c tpl.h
	gcc -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wa,-q tpl.c -o tpl.o

clean:
	rm tpl.o
	rm ../../lib/libtpl.a
