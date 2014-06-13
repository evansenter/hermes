include ../../Makefile.inc

HEADER   = h
CODE     = c
GCCFLAGS = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) $(GLOBAL_HEADERS)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	GCCFLAGS += -Ofast -march=native
else
	GCCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	GCCFLAGS += -fdiagnostics-color=always
endif

$(LIB)/libmulti_param.a: $(CODE)/multi_param.o
	ar cr $(LIB)/libmulti_param.a $(CODE)/multi_param.o

$(CODE)/multi_param.o: $(CODE)/multi_param.c $(HEADER)/multi_param.h
	$(GCC) $(GCCFLAGS) $(CODE)/multi_param.c -o $(CODE)/multi_param.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libmulti_param.a

install: $(LIB)/libmulti_param.a
	cp $(LIB)/libmulti_param.a $(LIBDIR)
