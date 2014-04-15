#ifndef SPECTRAL_FUNCTIONS_H
#define SPECTRAL_FUNCTIONS_H

#include "shared/libmfpt_header.h"
#include "data_structures.h"

EIGENSYSTEM eigensystem_from_row_ordered_transition_matrix(TRANSITION_MATRIX);
void population_proportion_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS, TRANSITION_MATRIX);
void equilibrium_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS, TRANSITION_MATRIX);
TRANSITION_MATRIX convert_structures_to_transition_matrix(const SOLUTION*, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(TRANSITION_MATRIX);
void invert_matrix(EIGENSYSTEM*);
double probability_at_time(const EIGENSYSTEM, double, int, int);
void find_key_structure_indices_in_structure_list(SPECTRAL_PARAMS*, const SOLUTION*, int, char*, char*);
void serialize_eigensystem(const EIGENSYSTEM, const SPECTRAL_PARAMS);
EIGENSYSTEM deserialize_eigensystem(const SPECTRAL_PARAMS);
long double estimate_equilibrium(const EIGENSYSTEM, const SPECTRAL_PARAMS);
int estimate_starting_index_to_scan_for_equilibrium(int, const EIGENSYSTEM, const SPECTRAL_PARAMS);
void print_equilibrium(const EIGENSYSTEM, const SPECTRAL_PARAMS);
void print_population_proportion(const EIGENSYSTEM, const SPECTRAL_PARAMS);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);
void print_eigenvalues(const EIGENSYSTEM);

#endif
