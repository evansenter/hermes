#ifndef FFTBOR2D_PARAMS_H
#define FFTBOR2D_PARAMS_H

#include "fftbor2d_data_structures.h"

FFTBOR2D_PARAMS init_fftbor2d_params();
void free_fftbor2d_params(FFTBOR2D_PARAMS);
void parse_fftbor2d_args(FFTBOR2D_PARAMS&, int, char**, void (*)(int));
void parse_fftbor2d_sequence_data(int, char**, int, FFTBOR2D_PARAMS&, void (*)(int));
char* find_energy_file(char*);
int fftbor2d_error_handling(const FFTBOR2D_PARAMS);
void debug_fftbor2d_parameters(const FFTBOR2D_PARAMS);
void fftbor2d_usage(int);

#endif
