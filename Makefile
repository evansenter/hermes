# Base Makefile

SUB_MAKE = make -C

all: tpl FFTbor2D RNAmfpt RNAeq MultiParam FFTmfpt FFTeq RateEq

install: all
	cp bin/FFTbor2D bin/RNAmfpt bin/RNAeq bin/FFTmfpt bin/FFTeq bin/RateEq ~/bin

tpl:
	$(SUB_MAKE) src/tpl

FFTbor2D:
	$(SUB_MAKE) src/fftbor2d
	cp src/fftbor2d/FFTbor2D.out bin/FFTbor2D

RNAmfpt:
	$(SUB_MAKE) src/mfpt
	cp src/mfpt/RNAmfpt.out bin/RNAmfpt

RNAeq:
	$(SUB_MAKE) src/population
	cp src/population/RNAeq.out bin/RNAeq

MultiParam:
	$(SUB_MAKE) src/multi_param

FFTmfpt:
	$(SUB_MAKE) mashup/mfpt_from_fftbor2d
	cp mashup/mfpt_from_fftbor2d/FFTmfpt.out bin/FFTmfpt

FFTeq:
	$(SUB_MAKE) mashup/population_from_fftbor2d
	cp mashup/population_from_fftbor2d/FFTeq.out bin/FFTeq

RateEq:
	$(SUB_MAKE) mashup/population_from_rate_matrix
	cp mashup/population_from_rate_matrix/RateEq.out bin/RateEq

clean:
	rm lib/libfftbor2d.a lib/libmfpt.a lib/libmulti_param.a lib/libpopulation.a lib/libtpl.a
	$(SUB_MAKE) src/tpl clean
	$(SUB_MAKE) src/fftbor2d clean
	$(SUB_MAKE) src/mfpt clean
	$(SUB_MAKE) src/population clean
	$(SUB_MAKE) src/multi_param clean
	$(SUB_MAKE) mashup/mfpt_from_fftbor2d clean
	$(SUB_MAKE) mashup/population_from_fftbor2d clean
	$(SUB_MAKE) mashup/population_from_rate_matrix clean
	rm bin/RNAmfpt bin/RNAeq bin/FFTbor2D bin/FFTmfpt bin/FFTeq bin/RateEq
