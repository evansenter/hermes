#ifndef KLP_MATRIX_PARSER_H
#define KLP_MATRIX_PARSER_H

#include "klp_matrix_data_structures.h"

KLP_MATRIX klp_matrix_from_file(const char*, short, void (*)());
TRANSITION_MATRIX transition_matrix_from_klp_matrix(KLP_MATRIX*, char);
int count_lines(const char*);
void populate_arrays(KLP_MATRIX*, const char*, short);

#endif
