#ifndef MFPT_ENERGY_GRID_H
#define MFPT_ENERGY_GRID_H

#include "lapack_externs.h"
#include "mfpt_data_structures.h"

TRANSITION_MATRIX convert_klp_matrix_to_transition_matrix(KLP_MATRIX*, MFPT_PARAMS*);
double compute_mfpt(const TRANSITION_MATRIX, const MFPT_PARAMS);
TRANSITION_MATRIX inverse(TRANSITION_MATRIX);
TRANSITION_MATRIX pseudoinverse(TRANSITION_MATRIX);
int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX*, MFPT_PARAMS*);
void set_bp_dist_from_start_and_end_positions(const KLP_MATRIX, MFPT_PARAMS*, int);
void extend_klp_matrix_to_all_possible_positions(KLP_MATRIX*, const MFPT_PARAMS);
void populate_remaining_probabilities_in_klp_matrix(KLP_MATRIX*, const MFPT_PARAMS);
double* populate_number_of_adjacent_moves(const KLP_MATRIX, const MFPT_PARAMS);
int number_of_permissible_single_bp_moves(const KLP_MATRIX, int);
TRANSITION_MATRIX populate_transition_matrix_from_stationary_matrix(const KLP_MATRIX, const MFPT_PARAMS, const double*, transition_probability);
double transition_rate_from_probabilities(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_probabilities_with_hastings(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies_with_hastings(const KLP_MATRIX, const double*, int, int, short);

#endif
