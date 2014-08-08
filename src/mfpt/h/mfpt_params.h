#ifndef MFPT_PARAMS_H
#define MFPT_PARAMS_H

#include "mfpt_data_structures.h"

MFPT_PARAMS init_mfpt_params();
void parse_mfpt_args(KLP_PARAMS*, MFPT_PARAMS*, int, char**, void (*)());
int mfpt_error_handling(const KLP_PARAMS, const MFPT_PARAMS);
void debug_mfpt_parameters(const MFPT_PARAMS);
void mfpt_usage();
void mfpt_flags();

#endif
