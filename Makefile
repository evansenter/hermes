# Base Makefile

all: RNAmfpt RNAspectral FFTbor2D MultiParam RNApopulation FFTmfpt RNAeq

install: all
	cp bin/RNAmfpt bin/RNAspectral bin/FFTbor2D bin/RNApopulation bin/FFTmfpt bin/RNAeq ~/bin

RNAmfpt:
	cd src/mfpt; make
	cp src/mfpt/RNAmfpt.out bin/RNAmfpt

RNAspectral:
	cd src/spectral; make
	cp src/spectral/RNAspectral.out bin/RNAspectral

FFTbor2D:
	cd src/fftbor2d; make
	cp src/fftbor2d/FFTbor2D.out bin/FFTbor2D

MultiParam:
	cd src/multi_param; make

RNApopulation:
	cd mashup/population_from_energy_grid; make
	cp mashup/population_from_energy_grid/RNApopulation.out bin/RNApopulation

FFTmfpt:
	cd mashup/mfpt_from_energy_grid; make
	cp mashup/mfpt_from_energy_grid/FFTmfpt.out bin/FFTmfpt

RNAeq:
	cd mashup/equilibrium_from_transition_rate_matrix; make
	cp mashup/equilibrium_from_transition_rate_matrix/RNAeq.out bin/RNAeq

clean:
	cd src/mfpt; make clean
	cd src/spectral; make clean
	cd src/fftbor2d; make clean
	cd src/multi_param; make clean
	cd mashup/population_from_energy_grid; make clean
	cd mashup/mfpt_from_energy_grid; make clean
	cd mashup/equilibrium_from_transition_rate_matrix; make clean
	rm bin/RNAmfpt bin/RNAspectral bin/FFTbor2D bin/RNApopulation bin/FFTmfpt bin/RNAeq
