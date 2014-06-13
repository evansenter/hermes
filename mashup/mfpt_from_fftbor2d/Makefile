include ../../Makefile.inc

GXXFLAGS = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) `pkg-config --cflags gsl`
LDFLAGS  = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack `pkg-config --libs gsl` -lfftbor2d -lmfpt -lmulti_param -lRNA

ifeq "$(GXX_GTEQ_4.6.0)" "1"
	GXXFLAGS += -Ofast -march=native
else
	GXXFLAGS += -O3
endif

ifeq "$(GXX_GTEQ_4.9.0)" "1"
	GXXFLAGS += -fdiagnostics-color=always
endif

FFTmfpt.out: $(CODE)/mfpt_from_fftbor2d.o
	$(MAKE) -C ../../src/fftbor2d
	$(MAKE) -C ../../src/mfpt
	$(GXX) -o FFTmfpt.out $(CODE)/mfpt_from_fftbor2d.o $(LDFLAGS)

$(CODE)/mfpt_from_fftbor2d.o: $(CODE)/mfpt_from_fftbor2d.cpp $(SHARED_HEADER)/shared/libfftbor2d_header.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(GXX) $(GXXFLAGS) $(CODE)/mfpt_from_fftbor2d.cpp -o $(CODE)/mfpt_from_fftbor2d.o

install: FFTmfpt.out
	cp FFTmfpt.out $(BINDIR)/FFTmfpt

clean:
	rm -f $(CODE)/*.o FFTmfpt.out
