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
	
RNAspectral.out: $(CODE)/rna_spectral.o $(CODE)/spectral_params.o $(CODE)/spectral_functions.o $(CODE)/spectral_initializers.o $(LIB)/libmfpt.a 
	$(CC) $(CODE)/rna_spectral.o $(CODE)/spectral_params.o $(CODE)/spectral_functions.o $(CODE)/spectral_initializers.o $(LDFLAGS) RNAspectral.out
	ar cr $(LIB)/libspectral.a $(CODE)/spectral_functions.o $(CODE)/spectral_params.o $(CODE)/spectral_initializers.o
	
$(CODE)/rna_spectral.o: $(CODE)/rna_spectral.c $(SHARED_HEADER)/vienna/data_structures.h $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/constants.h $(HEADER)/initializers.h
	$(CC) $(CCFLAGS) $(CODE)/rna_spectral.c -o $(CODE)/rna_spectral.o
	
$(CODE)/spectral_functions.o: $(CODE)/spectral_functions.c $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/data_structures.h $(HEADER)/initializers.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(CC) $(CCFLAGS) $(CODE)/spectral_functions.c -o $(CODE)/spectral_functions.o
	
$(CODE)/spectral_params.o: $(CODE)/spectral_params.c $(HEADER)/params.h $(HEADER)/data_structures.h
	$(CC) $(CCFLAGS) $(CODE)/spectral_params.c -o $(CODE)/spectral_params.o
  
$(CODE)/spectral_initializers.o: $(CODE)/spectral_initializers.c $(HEADER)/initializers.h $(HEADER)/functions.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(CC) $(CCFLAGS) $(CODE)/spectral_initializers.c -o $(CODE)/spectral_initializers.o
  
$(LIB)/libmfpt.a:
	cd ../mfpt; make

clean:
	rm -f $(CODE)/*.o $(LIB)/libspectral.a RNAspectral.out
	
install: RNAspectral.out
	cp RNAspectral.out $(BINDIR)/RNAspectral
	cp $(LIB)/libspectral.a $(LIBDIR)
		
