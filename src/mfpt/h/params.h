#ifndef MFPT_PARAMS_H
#define MFPT_PARAMS_H

#include "data_structures.h"

MFPT_PARAMS init_mfpt_params();
void parse_mfpt_args(MFPT_PARAMS*, int, char**);
int mfpt_error_handling(const MFPT_PARAMS);
void debug_mfpt_parameters(const MFPT_PARAMS);
void mfpt_usage();

#endif
