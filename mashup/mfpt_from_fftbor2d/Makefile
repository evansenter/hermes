include ../../Makefile.inc

GXXFLAGS = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER)
LDFLAGS  = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lfftbor2d -lmfpt -lmulti_param -lRNA -o

ifeq "$(GXX_GTEQ_4.6.0)" "1"
	GXXFLAGS += -Ofast -march=native
else
	GXXFLAGS += -O3
endif

ifeq "$(GXX_GTEQ_4.9.0)" "1"
	GXXFLAGS += -fdiagnostics-color=always
endif

FFTmfpt.out: $(CODE)/mfpt_from_fftbor2d.o
	cd ../../src/fftbor2d; make
	cd ../../src/mfpt; make
	$(GXX) $(CODE)/mfpt_from_fftbor2d.o $(LDFLAGS) FFTmfpt.out

$(CODE)/mfpt_from_fftbor2d.o: $(CODE)/mfpt_from_fftbor2d.cpp $(SHARED_HEADER)/shared/libfftbor2d_header.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(GXX) $(GXXFLAGS) $(CODE)/mfpt_from_fftbor2d.cpp -o $(CODE)/mfpt_from_fftbor2d.o

install: FFTmfpt.out
	cp FFTmfpt.out $(BINDIR)/FFTmfpt

clean:
	rm -f $(CODE)/*.o FFTmfpt.out
