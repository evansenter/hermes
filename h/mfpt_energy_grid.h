#ifndef MFPT_ENERGY_GRID_H
#define MFPT_ENERGY_GRID_H

#include "mfpt_data_structures.h"

double** convert_energy_grid_to_transition_matrix(int**, int**, double**, int*, MFPT_PARAMETERS);
double compute_mfpt(int*, int*, double**, int, MFPT_PARAMETERS);
double* inverse(double*, int);
double* pseudoinverse(double*, int);
int number_of_permissible_single_bp_moves(int, int, int*, int*, int);
double transition_rate_from_probabilities(double, double, double);
double transition_rate_from_energies(double, double, double);
double transition_rate_from_probabilities_with_hastings(double, double, double, double);
double transition_rate_from_energies_with_hastings(double, double, double, double);

#endif
