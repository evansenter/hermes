# Base Makefile

all: FFTbor2D CSVmfpt RNApopulation MultiParam mfpt_from_fftbor2d pop_from_fftbor2d pop_from_rate_matrix

install: all
	cp bin/FFTbor2D bin/CSVmfpt bin/RNApopulation bin/mfpt_from_fftbor2d bin/pop_from_fftbor2d bin/pop_from_rate_matrix ~/bin

FFTbor2D:
	cd src/fftbor2d; make
	cp src/fftbor2d/FFTbor2D.out bin/FFTbor2D

CSVmfpt:
	cd src/mfpt; make
	cp src/mfpt/CSVmfpt.out bin/CSVmfpt

RNApopulation:
	cd src/population; make
	cp src/population/RNApopulation.out bin/RNApopulation

MultiParam:
	cd src/multi_param; make

mfpt_from_fftbor2d:
	cd mashup/mfpt_from_fftbor2d; make
	cp mashup/mfpt_from_fftbor2d/mfpt_from_fftbor2d.out bin/mfpt_from_fftbor2d

pop_from_fftbor2d:
	cd mashup/population_from_fftbor2d; make
	cp mashup/population_from_fftbor2d/pop_from_fftbor2d.out bin/pop_from_fftbor2d

pop_from_rate_matrix:
	cd mashup/population_from_rate_matrix; make
	cp mashup/population_from_rate_matrix/pop_from_rate_matrix.out bin/pop_from_rate_matrix

clean:
	cd src/fftbor2d; make clean
	cd src/mfpt; make clean
	cd src/population; make clean
	cd src/multi_param; make clean
	cd mashup/mfpt_from_fftbor2d; make clean
	cd mashup/population_from_fftbor2d; make clean
	cd mashup/population_from_rate_matrix; make clean
	rm bin/CSVmfpt bin/RNApopulation bin/FFTbor2D bin/mfpt_from_fftbor2d bin/pop_from_fftbor2d bin/pop_from_rate_matrix
