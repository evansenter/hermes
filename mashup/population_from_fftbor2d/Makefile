include ../../Makefile.inc

GXXFLAGS = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER)
LDFLAGS  = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lfftbor2d -lmfpt -lpopulation -lmulti_param -ltpl -lRNA -o

ifeq "$(GXX_GTEQ_4.6.0)" "1"
	GXXFLAGS += -Ofast -march=native
else
	GXXFLAGS += -O3
endif

ifeq "$(GXX_GTEQ_4.9.0)" "1"
	GXXFLAGS += -fdiagnostics-color=always
endif

FFTeq.out: $(CODE)/population_from_fftbor2d.o
	cd ../../src/fftbor2d; make
	cd ../../src/mfpt; make
	cd ../../src/population; make
	$(GXX) $(CODE)/population_from_fftbor2d.o $(LDFLAGS) FFTeq.out

$(CODE)/population_from_fftbor2d.o: $(CODE)/population_from_fftbor2d.cpp $(SHARED_HEADER)/shared/libfftbor2d_header.h $(SHARED_HEADER)/shared/libmfpt_header.h $(SHARED_HEADER)/shared/libpopulation_header.h  $(SHARED_HEADER)/shared/libmulti_param_header.h
	$(GXX) $(GXXFLAGS) $(CODE)/population_from_fftbor2d.cpp -o $(CODE)/population_from_fftbor2d.o

install: FFTeq.out
	cp FFTeq.out $(BINDIR)/FFTeq

clean:
	rm -f $(CODE)/*.o FFTeq.out
