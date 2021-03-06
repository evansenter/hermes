#ifndef FFTBOR2D_PARAMS_H
#define FFTBOR2D_PARAMS_H

#include "fftbor2d_data_structures.h"

FFTBOR2D_PARAMS init_fftbor2d_params();
void free_fftbor2d_params(FFTBOR2D_PARAMS);
void parse_fftbor2d_args(FFTBOR2D_PARAMS&, int, char**, void (*)());
void parse_fftbor2d_data_from_file(FFTBOR2D_PARAMS&, void (*)());
void sanitize_sequence_and_structure_parameters(FFTBOR2D_PARAMS&);
int fftbor2d_error_handling(const FFTBOR2D_PARAMS);
void debug_fftbor2d_parameters(const FFTBOR2D_PARAMS);
void fftbor2d_usage();
void fftbor2d_flags();

#endif
