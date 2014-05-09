# Base Makefile

all: tpl FFTbor2D RNAmfpt RNAeq MultiParam FFTmfpt FFTeq RateEq

install: all
	cp bin/FFTbor2D bin/RNAmfpt bin/RNAeq bin/FFTmfpt bin/FFTeq bin/RateEq ~/bin

tpl:
	cd src/tpl; make

FFTbor2D:
	cd src/fftbor2d; make
	cp src/fftbor2d/FFTbor2D.out bin/FFTbor2D

RNAmfpt:
	cd src/mfpt; make
	cp src/mfpt/RNAmfpt.out bin/RNAmfpt

RNAeq:
	cd src/population; make
	cp src/population/RNAeq.out bin/RNAeq

MultiParam:
	cd src/multi_param; make

FFTmfpt:
	cd mashup/mfpt_from_fftbor2d; make
	cp mashup/mfpt_from_fftbor2d/FFTmfpt.out bin/FFTmfpt

FFTeq:
	cd mashup/population_from_fftbor2d; make
	cp mashup/population_from_fftbor2d/FFTeq.out bin/FFTeq

RateEq:
	cd mashup/population_from_rate_matrix; make
	cp mashup/population_from_rate_matrix/RateEq.out bin/RateEq

clean:
	rm -f lib/*
	cd src/tpl; make clean
	cd src/fftbor2d; make clean
	cd src/mfpt; make clean
	cd src/population; make clean
	cd src/multi_param; make clean
	cd mashup/mfpt_from_fftbor2d; make clean
	cd mashup/population_from_fftbor2d; make clean
	cd mashup/population_from_rate_matrix; make clean
	rm bin/RNAmfpt bin/RNAeq bin/FFTbor2D bin/FFTmfpt bin/FFTeq bin/RateEq
