hermes
======

All your kinetics are belong to us.
-----------------------------------

    Name                    Role
    --------                ----
    FFTbor2D                (2D energy landscape parameterized by bp. distance from two input structures for given sequence)
    RNAmfpt                 (Mean first passage time [hitting time] for input matrix in CSV (i, j, x) format)
    RNAeq                   (Population proportion and equilibrium for suboptimal structures of given RNA using spectral decomposition)
    FFTmfpt                 ([--fftbor2d --mfpt] Approximate mean first passage time [hitting time] using 2D energy landscape from FFTbor2D)
    FFTeq                   ([--fftbor2d --mfpt --population] Approximate population proportion and equilibrium using 2D energy landscape from FFTbor2D)
    RateEq                  ([--mfpt --population] Population proportion and equilibrium using row-ordered rate matrix in CSV (i, j, x) format)
    
Dependencies:
-------------
  
* CMake >= 2.6-patch 4, tested through 3.0.0
* GNU99 compiler, developed with gcc-4.9
* C++98 compiler supporting OpenMP, developed with g++-4.9
* GSL <= 1.15, tested through 1.16
* FFTW3 >= 3.3.3, tested through 3.3.4
* libRNA.a >= 2.0.7, tested through 2.1.7
  

Installation:
-------------

    cd build && cmake .. && make

Useful CMake configuration flags:
---------------------------------

    -DCMAKE_C_COMPILER=/path/to/c/compiler
    -DCMAKE_CXX_COMPILER=/path/to/c++/compiler
    -DCMAKE_INCLUDE_PATH=~/additional/headers
    -DCMAKE_LIBRARY_PATH="~/additional/libraries;/separated/by/semicolons;/and/quoted"
    -DCMAKE_INSTALL_PREFIX=/make/install/path/prefix

Note--this project uses the following commands to maintain homogenous styling:

    astyle --style=google --indent=spaces=2 --indent-switches --indent-cases --indent-namespaces --indent-labels --indent-col1-comments --break-blocks --pad-oper --pad-header --unpad-paren --fill-empty-lines --align-pointer=type --add-brackets --convert-tabs --recursive *.cpp *.c *.h
