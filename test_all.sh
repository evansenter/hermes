#!/usr/bin/env bash

echo "Testing FFTbor2D, expecting:
0	3	+0.10047371	+0.00000002
1	2	+0.00019384	+3.85251038
1	4	+0.00181492	+2.47390412
2	1	+0.02666435	+0.81761929
2	3	+0.00036362	+3.46479324
2	5	+0.16642614	-0.31104020
3	0	+0.70406341	-1.20000002"
echo "Result:"
FFTbor2D -s GGGAAACCC "........." "(((...)))"

echo "Testing RNAmfpt, expecting:
+699.65148561"
echo "Result:"
RNAmfpt -c ./src/mfpt/example.csv -xh

echo "Testing RNAeq, expecting:
+1.000000	+0.70407025	+0.10047166"
echo "Result:"
RNAeq -o -s GGGAAACCC -l "(((...)))" -i 0 -j 1 -p 1

echo "Testing FFTmfpt, expecting:
6786.738209"
echo "Result:"
FFTmfpt --fftbor2d-i GGGAAACCC --fftbor2d-j "........." --fftbor2d-k "(((...)))" --mfpt-x --mfpt-h

echo "Testing FFTeq, expecting:
+3.000000	+0.33271632	+0.03369467"
echo "Result:"
FFTeq --fftbor2d-i GGGGGCCCCC --fftbor2d-j ".........." --fftbor2d-k "(((....)))" --population-i 2 --population-j 3 --population-p 1

echo "Testing RateEq, expecting:
3.235000"
echo "Result:"
RateEq --mfpt-c ./mashup/population_from_rate_matrix/all_structures_for_gggaaaccc__19_5_rate_move_without_hastings.csv --population-a 19 --population-z 5 --population-e 1e-4
