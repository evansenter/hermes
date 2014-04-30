hermes
======

All your kinetics are belong to us.

Original name     New name                Role
-------------     --------                ----
FFTbor2D          FFTbor2D                (2D energy landscape parameterized by bp. distance from two input structures for given sequence)
RNAmfpt           CSVmfpt                 (Mean first passage time [hitting time] for input matrix in CSV (i, j, x) format)
RNAspectral       RNApopulation           (Population proportion and equilibrium for suboptimal structures of given RNA using spectral decomposition)
FFTmfpt           mfpt_from_fftbor2d      ([--fftbor2d --mfpt] Approximate mean first passage time [hitting time] using 2D energy landscape from FFTbor2D)
RNApopulation     pop_from_fftbor2d       ([--fftbor2d --mfpt --population] Approximate population proportion and equilibrium using 2D energy landscape from FFTbor2D)
RNAeq             pop_from_rate_matrix    ([--mfpt --population] Population proportion and equilibrium using row-ordered rate matrix in CSV (i, j, x) format)

Note--This project uses the following commands:

    astyle --style=google --indent=spaces=2 --indent-switches --indent-cases --indent-namespaces --indent-labels --indent-col1-comments --break-blocks --pad-oper --pad-header --unpad-paren --fill-empty-lines --align-pointer=type --add-brackets --convert-tabs --recursive *.cpp *.c *.h
