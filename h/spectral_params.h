#ifndef SPECTRAL_PARAMS_H
#define SPECTRAL_PARAMS_H

#include "spectral_data_structures.h"

SPECTRAL_PARAMS init_spectral_params();
SPECTRAL_PARAMS parse_spectral_args(int, char*[]);
int spectral_error_handling(SPECTRAL_PARAMS);
void debug_spectral_parameters(SPECTRAL_PARAMS);
void spectral_usage();

#endif
