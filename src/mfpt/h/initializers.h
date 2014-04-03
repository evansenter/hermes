#ifndef MFPT_INITIALIZERS_H
#define MFPT_INITIALIZERS_H

#include "data_structures.h"

#define ROW_ORDER(matrix, i, j, n) ((matrix)[((i) * (n)) + (j)])
#define COL_ORDER(matrix, i, j, n) ((matrix)[((j) * (n)) + (i)])

KLP_MATRIX init_klp_matrix(int);
void free_klp_matrix(KLP_MATRIX);
void print_klp_matrix(KLP_MATRIX);
TRANSITION_MATRIX init_transition_matrix(int, char);
TRANSITION_MATRIX transpose_matrix(TRANSITION_MATRIX);
void free_transition_matrix(TRANSITION_MATRIX);
void print_transition_matrix(KLP_MATRIX, double*, MFPT_PARAMS);

#endif
