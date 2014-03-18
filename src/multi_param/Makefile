# Makefile for RNAspectral

CCFLAGS          = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) -I ~/local/include
LDFLAGS          = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lm -lgomp -llapack -lgslcblas -lgsl -lmfpt -lRNA -o
BINDIR           = ~/bin
LIBDIR           = ~/lib
CC               = gcc
LIB              = ../../lib
SHARED_HEADER    = ../../h
HEADER           = h
CODE             = c
GCC_VERSION      = $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  = $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40900)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif

$(LIB)/libmulti_param.a: $(CODE)/multi_param.o
	ar cr $(LIB)/libmulti_param.a $(CODE)/multi_param.o

$(CODE)/multi_param.o: $(CODE)/multi_param.c $(HEADER)/multi_param.h
	$(CC) $(CCFLAGS) $(CODE)/multi_param.c -o $(CODE)/multi_param.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libmulti_param.a

install: $(LIB)/libmulti_param.a
	cp $(LIB)/libmulti_param.a $(LIBDIR)
