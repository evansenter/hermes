#ifndef KLP_MATRIX_FUNCTIONS_H
#define KLP_MATRIX_FUNCTIONS_H

#include "klp_matrix_data_structures.h"

typedef double(*transition_probability)(const KLP_MATRIX, const double*, int, int, short);

TRANSITION_MATRIX transpose_matrix(TRANSITION_MATRIX);
TRANSITION_MATRIX inverse(TRANSITION_MATRIX);
TRANSITION_MATRIX convert_klp_matrix_to_transition_matrix(KLP_MATRIX*, KLP_PARAMS*);
int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX*, KLP_PARAMS*);
void set_bp_dist_from_start_and_end_positions(const KLP_MATRIX, KLP_PARAMS*, int);
void extend_klp_matrix_to_all_possible_positions(KLP_MATRIX*, const KLP_PARAMS);
void populate_remaining_probabilities_in_klp_matrix(KLP_MATRIX*, const KLP_PARAMS);
double* populate_number_of_adjacent_moves(const KLP_MATRIX, const KLP_PARAMS);
int number_of_permissible_single_bp_moves(const KLP_MATRIX, int);
TRANSITION_MATRIX populate_transition_matrix_from_stationary_matrix(const KLP_MATRIX, const KLP_PARAMS, const double*, transition_probability);
double transition_rate_from_probabilities(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_probabilities_with_hastings(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies_with_hastings(const KLP_MATRIX, const double*, int, int, short);

#define ONE_BP_MOVE(i, j) ((int)abs(klp_matrix.k[(i)] - klp_matrix.k[(j)]) == 1 && (int)abs(klp_matrix.l[(i)] - klp_matrix.l[(j)]) == 1)
#define NONZERO_TO_NONZERO_PROB(i, j) (klp_matrix.p[(i)] > 0 && klp_matrix.p[(j)] > 0)

#endif