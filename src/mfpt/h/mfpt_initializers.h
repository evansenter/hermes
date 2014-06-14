#ifndef MFPT_INITIALIZERS_H
#define MFPT_INITIALIZERS_H

#include "mfpt_data_structures.h"

KLP_MATRIX init_klp_matrix(int);
void free_klp_matrix(KLP_MATRIX);
void print_klp_matrix(const KLP_MATRIX);
TRANSITION_MATRIX init_transition_matrix(int, char);
TRANSITION_MATRIX transpose_matrix(TRANSITION_MATRIX);
void free_transition_matrix(TRANSITION_MATRIX);
void print_transition_matrix(const KLP_MATRIX, const TRANSITION_MATRIX, const MFPT_PARAMS);

#endif
