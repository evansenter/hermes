# Makefile for RNAeq

CCFLAGS          = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) -I ~/local/include
LDFLAGS          = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lm -lgomp -llapack -lgslcblas -lgsl -ltpl -lmfpt -lRNA -o
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

RNAeq.out: $(CODE)/rna_population.o $(CODE)/population_params.o $(CODE)/population_functions.o $(CODE)/population_initializers.o
	cd ../mfpt; make
	$(CC) $(CODE)/rna_population.o $(CODE)/population_params.o $(CODE)/population_functions.o $(CODE)/population_initializers.o $(LDFLAGS) RNAeq.out
	ar cr $(LIB)/libpopulation.a $(CODE)/population_functions.o $(CODE)/population_params.o $(CODE)/population_initializers.o

$(CODE)/rna_population.o: $(CODE)/rna_population.c $(SHARED_HEADER)/vienna/data_structures.h $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/constants.h $(HEADER)/initializers.h
	$(CC) $(CCFLAGS) $(CODE)/rna_population.c -o $(CODE)/rna_population.o

$(CODE)/population_functions.o: $(CODE)/population_functions.c $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/data_structures.h $(HEADER)/initializers.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(CC) $(CCFLAGS) $(CODE)/population_functions.c -o $(CODE)/population_functions.o

$(CODE)/population_params.o: $(CODE)/population_params.c $(HEADER)/params.h $(HEADER)/data_structures.h
	$(CC) $(CCFLAGS) $(CODE)/population_params.c -o $(CODE)/population_params.o

$(CODE)/population_initializers.o: $(CODE)/population_initializers.c $(HEADER)/initializers.h $(HEADER)/functions.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(CC) $(CCFLAGS) $(CODE)/population_initializers.c -o $(CODE)/population_initializers.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libpopulation.a RNAeq.out

install: RNAeq.out
	cp RNAeq.out $(BINDIR)/RNAeq
	cp $(LIB)/libpopulation.a $(LIBDIR)
