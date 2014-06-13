include ../../Makefile.inc

GCCFLAGS = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) `pkg-config --cflags gsl`
LDFLAGS  = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack `pkg-config --libs gsl` -lmfpt -lpopulation -lmulti_param -ltpl -lRNA

ifeq "$(GGCC_GTEQ_4.6.0)" "1"
	GCCFLAGS += -Ofast -march=native
else
	GCCFLAGS += -O3
endif

ifeq "$(GGCC_GTEQ_4.9.0)" "1"
	GCCFLAGS += -fdiagnostics-color=always
endif

RateEq.out: $(CODE)/population_from_rate_matrix.o
	$(MAKE) -C ../../src/mfpt
	$(MAKE) -C ../../src/population
	$(GCC) -o RateEq.out $(CODE)/population_from_rate_matrix.o $(LDFLAGS)

$(CODE)/population_from_rate_matrix.o: $(CODE)/population_from_rate_matrix.c $(SHARED_HEADER)/shared/libmfpt_header.h $(SHARED_HEADER)/shared/libpopulation_header.h $(SHARED_HEADER)/shared/libmulti_param_header.h
	$(GCC) $(GCCFLAGS) $(CODE)/population_from_rate_matrix.c -o $(CODE)/population_from_rate_matrix.o

install: RateEq.out
	cp RateEq.out $(BINDIR)/RateEq

clean:
	rm -f $(CODE)/*.o RateEq.out
