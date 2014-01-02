# Base Makefile

all: src/mfpt/RNAmfpt src/spectral/RNAspectral src/fftbor2d/FFTbor2D
	
src/mfpt/RNAmfpt:
	cd src/mfpt; make
	
src/spectral/RNAspectral:
	cd src/spectral; make
	
src/fftbor2d/FFTbor2D:
	cd src/fftbor2d; make
	
clean:
	cd src/mfpt; make clean
	cd src/spectral; make clean
	cd src/fftbor2d; make clean
	
