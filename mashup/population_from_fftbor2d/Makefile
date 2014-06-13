include ../../Makefile.inc

GXXFLAGS = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) `pkg-config --cflags gsl`
LDFLAGS  = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack `pkg-config --libs gsl` -lfftbor2d -lmfpt -lpopulation -lmulti_param -ltpl -lRNA

ifeq "$(GXX_GTEQ_4.6.0)" "1"
	GXXFLAGS += -Ofast -march=native
else
	GXXFLAGS += -O3
endif

ifeq "$(GXX_GTEQ_4.9.0)" "1"
	GXXFLAGS += -fdiagnostics-color=always
endif

FFTeq.out: $(CODE)/population_from_fftbor2d.o
	$(MAKE) -C ../../src/fftbor2d
	$(MAKE) -C ../../src/mfpt
	$(MAKE) -C ../../src/population
	$(GXX) -o FFTeq.out $(CODE)/population_from_fftbor2d.o $(LDFLAGS)

$(CODE)/population_from_fftbor2d.o: $(CODE)/population_from_fftbor2d.cpp $(SHARED_HEADER)/shared/libfftbor2d_header.h $(SHARED_HEADER)/shared/libmfpt_header.h $(SHARED_HEADER)/shared/libpopulation_header.h  $(SHARED_HEADER)/shared/libmulti_param_header.h
	$(GXX) $(GXXFLAGS) $(CODE)/population_from_fftbor2d.cpp -o $(CODE)/population_from_fftbor2d.o

install: FFTeq.out
	cp FFTeq.out $(BINDIR)/FFTeq

clean:
	rm -f $(CODE)/*.o FFTeq.out
