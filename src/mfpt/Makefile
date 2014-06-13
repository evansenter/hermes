include ../../Makefile.inc

GCCFLAGS = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER) $(GLOBAL_HEADERS) `pkg-config --cflags gsl`
LDFLAGS  = -L . -L $(LIB) -L ~/lib $(GLOBAL_LIBS) -L /usr/local/include -lgomp -llapack `pkg-config --libs gsl`
HEADER   = h
CODE     = c

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	GCCFLAGS += -Ofast -march=native
else
	GCCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	GCCFLAGS += -fdiagnostics-color=always
endif

RNAmfpt.out: $(CODE)/csv_mfpt.o $(CODE)/mfpt_params.o $(CODE)/mfpt_functions.o $(CODE)/mfpt_parser.o $(CODE)/mfpt_initializers.o
	$(GCC) -o RNAmfpt.out $(CODE)/csv_mfpt.o $(CODE)/mfpt_parser.o $(CODE)/mfpt_params.o $(CODE)/mfpt_functions.o $(CODE)/mfpt_initializers.o $(LDFLAGS)
	ar cr $(LIB)/libmfpt.a $(CODE)/mfpt_functions.o $(CODE)/mfpt_params.o $(CODE)/mfpt_parser.o $(CODE)/mfpt_initializers.o

$(CODE)/csv_mfpt.o: $(CODE)/csv_mfpt.c $(HEADER)/parser.h $(HEADER)/functions.h $(HEADER)/constants.h $(HEADER)/params.h $(HEADER)/initializers.h
	$(GCC) $(GCCFLAGS) $(CODE)/csv_mfpt.c -o $(CODE)/csv_mfpt.o

$(CODE)/mfpt_initializers.o: $(CODE)/mfpt_initializers.c $(HEADER)/initializers.h
	$(GCC) $(GCCFLAGS) $(CODE)/mfpt_initializers.c -o $(CODE)/mfpt_initializers.o

$(CODE)/mfpt_parser.o: $(CODE)/mfpt_parser.c $(HEADER)/parser.h
	$(GCC) $(GCCFLAGS) $(CODE)/mfpt_parser.c -o $(CODE)/mfpt_parser.o

$(CODE)/mfpt_params.o: $(CODE)/mfpt_params.c $(HEADER)/params.h
	$(GCC) $(GCCFLAGS) $(CODE)/mfpt_params.c -o $(CODE)/mfpt_params.o

$(CODE)/mfpt_functions.o: $(CODE)/mfpt_functions.c $(HEADER)/functions.h $(HEADER)/constants.h $(HEADER)/params.h $(HEADER)/initializers.h
	$(GCC) $(GCCFLAGS) $(CODE)/mfpt_functions.c -o $(CODE)/mfpt_functions.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libmfpt.a RNAmfpt.out

install: RNAmfpt.out
	cp RNAmfpt $(BINDIR)/RNAmfpt
	cp $(LIB)/libmfpt.a $(LIBDIR)
