#ifndef KLP_MATRIX_PARAMS_H
#define KLP_MATRIX_PARAMS_H

#include "klp_matrix_data_structures.h"

KLP_PARAMS init_klp_matrix_params();
void parse_klp_matrix_args(KLP_PARAMS*, int, char**, void (*)());
int klp_matrix_error_handling(const KLP_PARAMS);
void debug_klp_matrix_parameters(const KLP_PARAMS);
void klp_matrix_usage();
void klp_matrix_flags();

#endif
