#ifndef MFPT_PARSER_H
#define MFPT_PARSER_H

#include "data_structures.h"

KLP_MATRIX klp_matrix_from_file(const MFPT_PARAMS);
double* transition_matrix_from_klp_matrix(KLP_MATRIX*);
int count_lines(char*);
void populate_arrays(KLP_MATRIX*, const MFPT_PARAMS);

#endif
