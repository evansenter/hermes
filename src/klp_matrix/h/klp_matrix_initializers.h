#ifndef KLP_MATRIX_INITIALIZERS_H
#define KLP_MATRIX_INITIALIZERS_H

#include "klp_matrix_data_structures.h"

KLP_MATRIX init_klp_matrix(int);
void free_klp_matrix(KLP_MATRIX);
void print_klp_matrix(const KLP_MATRIX);
TRANSITION_MATRIX init_transition_matrix(int, char);
void free_transition_matrix(TRANSITION_MATRIX);
void print_transition_matrix(const TRANSITION_MATRIX);
void print_transition_matrix_with_klp_positions(const KLP_MATRIX, const TRANSITION_MATRIX);

#endif
