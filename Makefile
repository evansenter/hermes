# Base Makefile

all: RNAmfpt RNAspectral FFTbor2D RNApopulation FFTmfpt

install: all
	cp RNAmfpt RNAspectral FFTbor2D RNApopulation FFTmfpt ~/bin
	
RNAmfpt:
	cd src/mfpt; make
	cp src/mfpt/RNAmfpt.out RNAmfpt
	
RNAspectral:
	cd src/spectral; make
	cp src/spectral/RNAspectral.out RNAspectral
	
FFTbor2D:
	cd src/fftbor2d; make
	cp src/fftbor2d/FFTbor2D.out FFTbor2D

RNApopulation:
	cd mashup/population_from_energy_grid; make
	cp mashup/population_from_energy_grid/RNApopulation.out RNApopulation
	
FFTmfpt:
	cd mashup/mfpt_from_energy_grid; make
	cp mashup/mfpt_from_energy_grid/FFTmfpt.out FFTmfpt
	
clean:
	cd src/mfpt; make clean
	cd src/spectral; make clean
	cd src/fftbor2d; make clean
	cd mashup/population_from_energy_grid; make clean
	cd mashup/mfpt_from_energy_grid; make clean
	rm RNAmfpt RNAspectral FFTbor2D RNApopulation FFTmfpt
