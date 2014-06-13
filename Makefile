include Makefile.inc

all: tpl FFTbor2D RNAmfpt RNAeq MultiParam FFTmfpt FFTeq RateEq

install: all
	cp bin/FFTbor2D bin/RNAmfpt bin/RNAeq bin/FFTmfpt bin/FFTeq bin/RateEq ~/bin

tpl:
	$(MAKE) -C src/tpl

FFTbor2D:
	$(MAKE) -C src/fftbor2d
	cp src/fftbor2d/FFTbor2D.out bin/FFTbor2D

RNAmfpt:
	$(MAKE) -C src/mfpt
	cp src/mfpt/RNAmfpt.out bin/RNAmfpt

RNAeq:
	$(MAKE) -C src/population
	cp src/population/RNAeq.out bin/RNAeq

MultiParam:
	$(MAKE) -C src/multi_param

FFTmfpt:
	$(MAKE) -C mashup/mfpt_from_fftbor2d
	cp mashup/mfpt_from_fftbor2d/FFTmfpt.out bin/FFTmfpt

FFTeq:
	$(MAKE) -C mashup/population_from_fftbor2d
	cp mashup/population_from_fftbor2d/FFTeq.out bin/FFTeq

RateEq:
	$(MAKE) -C mashup/population_from_rate_matrix
	cp mashup/population_from_rate_matrix/RateEq.out bin/RateEq

clean:
	rm -f lib/libfftbor2d.a lib/libmfpt.a lib/libmulti_param.a lib/libpopulation.a lib/libtpl.a
	$(MAKE) -C src/tpl clean
	$(MAKE) -C src/fftbor2d clean
	$(MAKE) -C src/mfpt clean
	$(MAKE) -C src/population clean
	$(MAKE) -C src/multi_param clean
	$(MAKE) -C mashup/mfpt_from_fftbor2d clean
	$(MAKE) -C mashup/population_from_fftbor2d clean
	$(MAKE) -C mashup/population_from_rate_matrix clean
	rm -f bin/RNAmfpt bin/RNAeq bin/FFTbor2D bin/FFTmfpt bin/FFTeq bin/RateEq
