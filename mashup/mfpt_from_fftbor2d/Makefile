CCFLAGS          = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wno-write-strings -Wa,-q -I $(HEADER) -I $(SHARED_HEADER)
LDFLAGS          = -L . -L $(LIB) -L ~/lib -L ~/local/lib -L /usr/local/include -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lfftbor2d -lmfpt -lmulti_param -lRNA -o
BINDIR           = ~/bin
LIBDIR           = ~/lib
CC               = g++
LIB              = ../../lib
SHARED_HEADER    = ../../h
HEADER           = .
CODE             = .
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

mfpt_from_fftbor2d.out: $(CODE)/mfpt_from_fftbor2d.o
	$(CC) $(CODE)/mfpt_from_fftbor2d.o $(LDFLAGS) mfpt_from_fftbor2d.out

$(CODE)/mfpt_from_fftbor2d.o: $(CODE)/mfpt_from_fftbor2d.cpp $(SHARED_HEADER)/shared/libfftbor2d_header.h $(SHARED_HEADER)/shared/libmfpt_header.h
	$(CC) $(CCFLAGS) $(CODE)/mfpt_from_fftbor2d.cpp -o $(CODE)/mfpt_from_fftbor2d.o

install: mfpt_from_fftbor2d.out
	cp mfpt_from_fftbor2d.out $(BINDIR)/mfpt_from_fftbor2d

clean:
	rm -f $(CODE)/*.o mfpt_from_fftbor2d.out
