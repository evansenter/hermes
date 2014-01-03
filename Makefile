# Base Makefile

all: src/mfpt/RNAmfpt src/spectral/RNAspectral src/fftbor2d/FFTbor2D mashup/population_from_energy_grid/RNApopulation mashup/mfpt_from_energy_grid/FFTmfpt

install: all
	cp RNAmfpt RNAspectral FFTbor2D RNApopulation FFTmfpt ~/bin
	
src/mfpt/RNAmfpt:
	cd src/mfpt; make
	cp src/mfpt/RNAmfpt.out RNAmfpt
	
src/spectral/RNAspectral:
	cd src/spectral; make
	cp src/spectral/RNAspectral.out RNAspectral
	
src/fftbor2d/FFTbor2D:
	cd src/fftbor2d; make
	cp src/fftbor2d/FFTbor2D.out FFTbor2D

mashup/population_from_energy_grid/RNApopulation:
	cd mashup/population_from_energy_grid; make
	cp mashup/population_from_energy_grid/RNApopulation.out RNApopulation
	
mashup/mfpt_from_energy_grid/FFTmfpt:
	cd mashup/mfpt_from_energy_grid; make
	cp mashup/mfpt_from_energy_grid/FFTmfpt.out FFTmfpt
	
clean:
	cd src/mfpt; make clean
	cd src/spectral; make clean
	cd src/fftbor2d; make clean
	cd mashup/population_from_energy_grid; make clean
	cd mashup/mfpt_from_energy_grid; make clean
	rm RNAmfpt RNAspectral FFTbor2D RNApopulation FFTmfpt
