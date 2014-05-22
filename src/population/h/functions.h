#ifndef POPULATION_FUNCTIONS_H
#define POPULATION_FUNCTIONS_H

#include "shared/libmfpt_header.h"
#include "data_structures.h"

void population_from_row_ordered_transition_matrix(POPULATION_PARAMS*, TRANSITION_MATRIX);
EIGENSYSTEM eigensystem_from_row_ordered_transition_matrix(TRANSITION_MATRIX);
void population_proportion_from_row_ordered_transition_matrix(const POPULATION_PARAMS, TRANSITION_MATRIX);
void equilibrium_from_row_ordered_transition_matrix(POPULATION_PARAMS*, TRANSITION_MATRIX);
TRANSITION_MATRIX convert_structures_to_transition_matrix(const SOLUTION*, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(TRANSITION_MATRIX);
void invert_matrix(EIGENSYSTEM*);
double probability_at_logtime(const EIGENSYSTEM, const POPULATION_PARAMS, double, int, int);
void find_key_structure_indices_in_structure_list(POPULATION_PARAMS*, const SOLUTION*, int);
void serialize_eigensystem(const EIGENSYSTEM, const POPULATION_PARAMS);
EIGENSYSTEM deserialize_eigensystem(const POPULATION_PARAMS);
void ensure_key_structures_and_energies_assigned(POPULATION_PARAMS*);
void set_mfe_str_and_energy(POPULATION_PARAMS*);
void set_energy_of_target_str(POPULATION_PARAMS*);
SOLUTION* sample_structures(const POPULATION_PARAMS);
double boltzmann_probability(const POPULATION_PARAMS);
double gsl_equilibrium_root_solver(double logtime, void* params);
double estimate_equilibrium(const EIGENSYSTEM, const POPULATION_PARAMS);
double soft_bound_for_population_proportion(const EIGENSYSTEM, const POPULATION_PARAMS, int, double, int);
int estimate_starting_index_to_scan_for_equilibrium(int, const EIGENSYSTEM, const POPULATION_PARAMS);
int index_in_equilibrium_within_window(const EIGENSYSTEM, const POPULATION_PARAMS, int, double);
void print_equilibrium(const EIGENSYSTEM, POPULATION_PARAMS*);
void print_population_proportion(const EIGENSYSTEM, const POPULATION_PARAMS);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);
void print_eigenvalues(const EIGENSYSTEM);

#endif
