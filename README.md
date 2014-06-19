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

Note--this project uses the following commands:

    astyle --style=google --indent=spaces=2 --indent-switches --indent-cases --indent-namespaces --indent-labels --indent-col1-comments --break-blocks --pad-oper --pad-header --unpad-paren --fill-empty-lines --align-pointer=type --add-brackets --convert-tabs --recursive *.cpp *.c *.h

Installation:
    cd build && cmake .. && make

Useful CMake flags:
    -DCMAKE_INSTALL_PREFIX:PATH=/install/path/prefix
    -DCMAKE_C_COMPILER=/path/to/c/compiler
    -DCMAKE_CXX_COMPILER=/path/to/c++/compiler
    -DCMAKE_INCLUDE_PATH:PATH=~/additional/headers
    -DCMAKE_LIBRARY_PATH:PATH="~/additional/libraries;/separated/by/semicolons;/and/quoted"

    (cd .. && git pull) && rm -rf * && cmake .. -DCMAKE_C_COMPILER:PATH=/cluster/home/evansenter/bin/gcc-4.8.1/bin/gcc -DCMAKE_CXX_COMPILER:PATH=/cluster/home/evansenter/bin/gcc-4.8.1/bin/g++ -DCMAKE_LIBRARY_PATH:PATH="/cluster/home/evansenter/lib;/cluster/home/evansenter/local/lib" -DCMAKE_INCLUDE_PATH:PATH=/cluster/home/evansenter/local/include && make