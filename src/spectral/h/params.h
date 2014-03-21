#ifndef SPECTRAL_PARAMS_H
#define SPECTRAL_PARAMS_H

#include "data_structures.h"

SPECTRAL_PARAMS init_spectral_params();
void parse_spectral_args(SPECTRAL_PARAMS*, int, char**);
int spectral_error_handling(const SPECTRAL_PARAMS);
void debug_spectral_parameters(const SPECTRAL_PARAMS);
void spectral_usage();

#endif
