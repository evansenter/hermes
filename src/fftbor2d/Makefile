include ../../Makefile.inc

GXXFLAGS = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) $(GLOBAL_HEADERS) `pkg-config --cflags gsl`
LDFLAGS  = -L . -L $(LIB) -L ~/lib $(GLOBAL_LIBS) -L /usr/local/include -lfftw3 -lgomp -llapack `pkg-config --libs gsl` -lRNA
HEADER   = h
CODE     = cpp

ifeq "$(GXX_GTEQ_4.6.0)" "1"
	GXXFLAGS += -Ofast -march=native
else
	GXXFLAGS += -O3
endif

ifeq "$(GXX_GTEQ_4.9.0)" "1"
	GXXFLAGS += -fdiagnostics-color=always
endif

FFTbor2D.out: $(CODE)/fftbor2d_functions.o $(CODE)/fftbor2d_initializers.o $(CODE)/fftbor2d_params.o $(CODE)/fftbor2d.o
	$(GXX) -o FFTbor2D.out $(CODE)/fftbor2d_functions.o $(CODE)/fftbor2d_initializers.o $(CODE)/fftbor2d_params.o $(CODE)/fftbor2d.o $(LDFLAGS)
	ar cr $(LIB)/libfftbor2d.a $(CODE)/fftbor2d_functions.o $(CODE)/fftbor2d_params.o $(CODE)/fftbor2d_initializers.o

$(CODE)/fftbor2d.o: $(CODE)/fftbor2d.cpp $(HEADER)/functions.h $(HEADER)/params.h $(HEADER)/initializers.h
	$(GXX) $(GXXFLAGS) $(CODE)/fftbor2d.cpp -o $(CODE)/fftbor2d.o

$(CODE)/fftbor2d_params.o: $(CODE)/fftbor2d_params.cpp $(HEADER)/params.h
	$(GXX) $(GXXFLAGS) $(CODE)/fftbor2d_params.cpp -o $(CODE)/fftbor2d_params.o

$(CODE)/fftbor2d_initializers.o: $(CODE)/fftbor2d_initializers.cpp $(HEADER)/initializers.h $(HEADER)/functions.h
	$(GXX) $(GXXFLAGS) $(CODE)/fftbor2d_initializers.cpp -o $(CODE)/fftbor2d_initializers.o

$(CODE)/fftbor2d_functions.o: $(CODE)/fftbor2d_functions.cpp $(HEADER)/functions.h $(HEADER)/initializers.h $(HEADER)/params.h
	$(GXX) $(GXXFLAGS) $(CODE)/fftbor2d_functions.cpp -o $(CODE)/fftbor2d_functions.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libfftbor2d.a FFTbor2D.out

install: FFTbor2D.out
	cp FFTbor2D.out $(BINDIR)/FFTbor2D
	cp $(LIB)/libfftbor2d.a $(LIBDIR)
	cp ../../misc/rna_turner1999.par ../../misc/rna_turner2004.par $(BINDIR)

