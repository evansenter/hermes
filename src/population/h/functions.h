#ifndef POPULATION_FUNCTIONS_H
#define POPULATION_FUNCTIONS_H

#include "shared/libmfpt_header.h"
#include "data_structures.h"

EIGENSYSTEM eigensystem_from_row_ordered_transition_matrix(TRANSITION_MATRIX);
void population_proportion_from_row_ordered_transition_matrix(const POPULATION_PARAMS, TRANSITION_MATRIX);
void equilibrium_from_row_ordered_transition_matrix(const POPULATION_PARAMS, TRANSITION_MATRIX);
TRANSITION_MATRIX convert_structures_to_transition_matrix(const SOLUTION*, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(TRANSITION_MATRIX);
void invert_matrix(EIGENSYSTEM*);
double probability_at_time(const EIGENSYSTEM, double, int, int);
void find_key_structure_indices_in_structure_list(POPULATION_PARAMS*, const SOLUTION*, int, char*, char*);
void serialize_eigensystem(const EIGENSYSTEM, const POPULATION_PARAMS);
EIGENSYSTEM deserialize_eigensystem(const POPULATION_PARAMS);
long double estimate_equilibrium(const EIGENSYSTEM, const POPULATION_PARAMS);
int estimate_starting_index_to_scan_for_equilibrium(int, const EIGENSYSTEM, const POPULATION_PARAMS);
int is_index_in_equilibrium_within_window_position(const EIGENSYSTEM, const POPULATION_PARAMS, int, int);
void print_equilibrium(const EIGENSYSTEM, const POPULATION_PARAMS);
void print_population_proportion(const EIGENSYSTEM, const POPULATION_PARAMS);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);
void print_eigenvalues(const EIGENSYSTEM);

#endif
