include ../../Makefile.inc

HEADER   = h
CODE     = c
GCCFLAGS = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) $(GLOBAL_HEADERS) `gsl-config --cflags`
LDFLAGS  = -L . -L $(LIB) $(GLOBAL_LIBS) -lgomp -llapack `gsl-config --libs` -ltpl -lmfpt -lRNA

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	GCCFLAGS += -Ofast -march=native
else
	GCCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	GCCFLAGS += -fdiagnostics-color=always
endif

RNAeq.out: $(CODE)/rna_population.o $(CODE)/population_params.o $(CODE)/population_functions.o $(CODE)/population_initializers.o
	$(MAKE) -C ../mfpt
	$(GCC) -o RNAeq.out $(CODE)/rna_population.o $(CODE)/population_params.o $(CODE)/population_functions.o $(CODE)/population_initializers.o $(LDFLAGS)
	ar cr $(LIB)/libpopulation.a $(CODE)/population_functions.o $(CODE)/population_params.o $(CODE)/population_initializers.o

$(CODE)/rna_population.o: $(CODE)/rna_population.c $(SHARED_HEADER)/vienna/data_structures.h $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/constants.h $(HEADER)/initializers.h
	$(GCC) $(GCCFLAGS) $(CODE)/rna_population.c -o $(CODE)/rna_population.o

$(CODE)/population_functions.o: $(CODE)/population_functions.c $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/data_structures.h $(HEADER)/initializers.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(GCC) $(GCCFLAGS) $(CODE)/population_functions.c -o $(CODE)/population_functions.o

$(CODE)/population_params.o: $(CODE)/population_params.c $(HEADER)/params.h $(HEADER)/data_structures.h
	$(GCC) $(GCCFLAGS) $(CODE)/population_params.c -o $(CODE)/population_params.o

$(CODE)/population_initializers.o: $(CODE)/population_initializers.c $(HEADER)/initializers.h $(HEADER)/functions.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(GCC) $(GCCFLAGS) $(CODE)/population_initializers.c -o $(CODE)/population_initializers.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libpopulation.a RNAeq.out

install: RNAeq.out
	cp RNAeq.out $(BINDIR)/RNAeq
	cp $(LIB)/libpopulation.a $(LIBDIR)
