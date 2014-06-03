#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "vienna/functions.h"
#include "shared/libmfpt_header.h"
#include "shared/libtpl_header.h"
#include "constants.h"
#include "params.h"
#include "initializers.h"
#include "functions.h"

void population_from_row_ordered_transition_matrix(const POPULATION_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
  if (parameters.equilibrium) {
    equilibrium_from_row_ordered_transition_matrix(parameters, row_transition_matrix);
  } else {
    population_proportion_from_row_ordered_transition_matrix(parameters, row_transition_matrix);
  }
}

EIGENSYSTEM eigensystem_from_row_ordered_transition_matrix(TRANSITION_MATRIX row_transition_matrix) {
  TRANSITION_MATRIX column_transition_matrix;
  EIGENSYSTEM eigensystem;
  
  column_transition_matrix = transpose_matrix(row_transition_matrix);
  eigensystem              = convert_transition_matrix_to_eigenvectors(column_transition_matrix);
  invert_matrix(&eigensystem);
  
  return eigensystem;
}

void population_proportion_from_row_ordered_transition_matrix(const POPULATION_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
  EIGENSYSTEM eigensystem;
  
  eigensystem = eigensystem_from_row_ordered_transition_matrix(row_transition_matrix);
  
  print_population_proportion(eigensystem, parameters);
}

void equilibrium_from_row_ordered_transition_matrix(POPULATION_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
  EIGENSYSTEM eigensystem;
  
  eigensystem = eigensystem_from_row_ordered_transition_matrix(row_transition_matrix);
  print_equilibrium(eigensystem, parameters);
  
}

TRANSITION_MATRIX convert_structures_to_transition_matrix(const SOLUTION* all_structures, int num_structures) {
  // This code to convert the list of structures to a transition matrix does not ensure that there are only single b.p. moves.
  int i, j;
  double col_sum;
  TRANSITION_MATRIX transition_matrix = init_transition_matrix(num_structures, 'R');
  
  for (i = 0; i < num_structures; ++i) {
    col_sum = 0;
    
    for (j = 0; j < num_structures; ++j) {
      if (i != j) {
        T_COL_ORDER(transition_matrix, i, j) = \
                                               MIN(1, exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT));
                                               
        col_sum += T_COL_ORDER(transition_matrix, i, j);
#ifdef INSANE_DEBUG
        printf("%d\t%d\t%.4e\n", i, j, T_COL_ORDER(transition_matrix, i, j));
#endif
      }
    }
    
#ifdef INSANE_DEBUG
    printf("%d col_sum:\t%.4e\n\n", i, col_sum);
#endif
    T_COL_ORDER(transition_matrix, i, i) = -col_sum;
  }
  
  return transition_matrix;
}

EIGENSYSTEM convert_transition_matrix_to_eigenvectors(TRANSITION_MATRIX transition_matrix) {
  // Generates right eigenvectors from a column-ordered transition rate matrix, and saves them in a one-dimensional column-ordered array (eigensystem.vectors).
  int i, j;
  EIGENSYSTEM eigensystem;
  eigensystem                             = init_eigensystem(transition_matrix.row_length);
  gsl_matrix_view matrix_view             = gsl_matrix_view_array(transition_matrix.matrix, eigensystem.length, eigensystem.length);
  gsl_vector_complex* eigenvalues         = gsl_vector_complex_alloc(eigensystem.length);
  gsl_matrix_complex* eigenvectors        = gsl_matrix_complex_alloc(eigensystem.length, eigensystem.length);
  gsl_eigen_nonsymmv_workspace* workspace = gsl_eigen_nonsymmv_alloc(eigensystem.length);
  gsl_eigen_nonsymmv_params(1, workspace);
  gsl_eigen_nonsymmv(&matrix_view.matrix, eigenvalues, eigenvectors, workspace);
  gsl_eigen_nonsymmv_free(workspace);
  
  for (i = 0; i < eigensystem.length; ++i) {
    eigensystem.values[i]               = GSL_REAL(gsl_vector_complex_get(eigenvalues, i));
    gsl_vector_complex_view eigenvector = gsl_matrix_complex_column(eigenvectors, i);
    
    for (j = 0; j < eigensystem.length; ++j) {
      E_COL_ORDER(eigensystem.vectors, i, j, eigensystem.length) = GSL_REAL(gsl_vector_complex_get(&eigenvector.vector, j));
    }
  }
  
  free_transition_matrix(transition_matrix);
  gsl_vector_complex_free(eigenvalues);
  gsl_matrix_complex_free(eigenvectors);
  return eigensystem;
}

void invert_matrix(EIGENSYSTEM* eigensystem) {
  int i, j, signum;
  gsl_matrix* matrix_to_invert = gsl_matrix_alloc(eigensystem->length, eigensystem->length);
  gsl_matrix* inversion_matrix = gsl_matrix_alloc(eigensystem->length, eigensystem->length);
  gsl_permutation* permutation = gsl_permutation_alloc(eigensystem->length);
  
  for (i = 0; i < eigensystem->length; ++i) {
    for (j = 0; j < eigensystem->length; ++j) {
      gsl_matrix_set(matrix_to_invert, i, j, E_ROW_ORDER(eigensystem->vectors, i, j, eigensystem->length));
    }
  }
  
  gsl_linalg_LU_decomp(matrix_to_invert, permutation, &signum);
  gsl_linalg_LU_invert(matrix_to_invert, permutation, inversion_matrix);
  
  for (i = 0; i < eigensystem->length; ++i) {
    for (j = 0; j < eigensystem->length; ++j) {
      E_ROW_ORDER(eigensystem->inverse_vectors, i, j, eigensystem->length) = gsl_matrix_get(inversion_matrix, i, j);
    }
  }
  
  gsl_matrix_free(matrix_to_invert);
  gsl_matrix_free(inversion_matrix);
  gsl_permutation_free(permutation);
}

double probability_at_logtime(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters, double timepoint, int start_index, int target_index) {
  // This function is hard-wired to only consider the kinetics for folding from a distribution where p_{0}(start_index) == 1.
  int i;
  double cumulative_probability = 0;
  
  if (parameters.start_index < 0 || parameters.end_index < 0) {
    if (parameters.start_index < 0) {
      fprintf(stderr, "Error: the starting structure could not be found. Usually this means that you did not specify an explicit starting structure so the empty structure was used, but the energy band was not wide enough for RNAsubopt to sample it.\n");
    }
    
    if (parameters.end_index < 0) {
      fprintf(stderr, "Error: the ending structure could not be found. Usually this means that you did explicitly provided a suboptimal ending structure, but the energy band was not wide enough for RNAsubopt to sample it.\n");
    }
    
    abort();
  }
  
  for (i = 0; i < eigensystem.length; ++i) {
    cumulative_probability +=
      E_ROW_ORDER(eigensystem.vectors, target_index, i, eigensystem.length) *
      E_ROW_ORDER(eigensystem.inverse_vectors, i, start_index, eigensystem.length) *
      exp(eigensystem.values[i] * pow(10, timepoint));
  }
  
  return cumulative_probability;
}

void find_key_structure_indices_in_structure_list(POPULATION_PARAMS* parameters, const SOLUTION* all_structures, int num_structures) {
  int i;
  
  for (i = 0; i < num_structures; ++i) {
    if (parameters->start_index == -1) {
      if (!strcmp(parameters->start_structure, all_structures[i].structure)) {
        parameters->start_index = i;
      }
    }
    
    if (parameters->end_index == -1) {
      if (!strcmp(parameters->end_structure, all_structures[i].structure)) {
        parameters->end_index = i;
      }
    }
  }
}

void serialize_eigensystem(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  tpl_node* tpl;
  int i;
  double value, vector, inverse_vector;
  tpl = tpl_map("iA(f)A(f)A(f)", &eigensystem.length, &value, &vector, &inverse_vector);
  tpl_pack(tpl, 0);
  
  for (i = 0; i < eigensystem.length * eigensystem.length; ++i) {
    if (i < eigensystem.length) {
      value = eigensystem.values[i];
      tpl_pack(tpl, 1);
    }
    
    vector         = eigensystem.vectors[i];
    inverse_vector = eigensystem.inverse_vectors[i];
    
    tpl_pack(tpl, 2);
    tpl_pack(tpl, 3);
  }
  
  tpl_dump(tpl, TPL_FILE, parameters.filename);
  tpl_free(tpl);
}

EIGENSYSTEM deserialize_eigensystem(const POPULATION_PARAMS parameters) {
  tpl_node* tpl;
  int i, length;
  double value, vector, inverse_vector;
  EIGENSYSTEM eigensystem;
  
  tpl = tpl_map("iA(f)A(f)A(f)", &length, &value, &vector, &inverse_vector);
  tpl_load(tpl, TPL_FILE, parameters.filename);
  tpl_unpack(tpl, 0);
  
  eigensystem = init_eigensystem(length);
  
  for (i = 0; i < eigensystem.length * eigensystem.length; ++i) {
    if (i < eigensystem.length) {
      tpl_unpack(tpl, 1);
      eigensystem.values[i] = value;
    }
    
    tpl_unpack(tpl, 2);
    tpl_unpack(tpl, 3);
    
    eigensystem.vectors[i]         = vector;
    eigensystem.inverse_vectors[i] = inverse_vector;
  }
  
  tpl_free(tpl);
  
  return eigensystem;
}

void ensure_key_structures_and_energies_assigned(POPULATION_PARAMS* parameters) {
  int i, seq_length = strlen(parameters->sequence);
  
  if (parameters->start_structure == NULL) {
    parameters->start_structure = malloc((seq_length + 1) * sizeof(char));
    
    for (i = 0; i < seq_length; ++i) {
      parameters->start_structure[i] = '.';
    }
    
    parameters->start_structure[i] = '\0';
  }
  
  if (parameters->end_structure == NULL) {
    set_end_structure(parameters);
  }
}

void set_end_structure(POPULATION_PARAMS* parameters) {
  paramT* vienna_params;
  
  vienna_params             = init_vienna_params(*parameters);
  parameters->end_structure = malloc((strlen(parameters->sequence) + 1) * sizeof(char));
  fold_par(parameters->sequence, parameters->end_structure, vienna_params, 0, 0);
}

SOLUTION* sample_structures(const POPULATION_PARAMS parameters) {
  double energy_cap;
  paramT* vienna_params;
  
  vienna_params = init_vienna_params(parameters);
  energy_cap    = parameters.energy_cap ? (int)round(parameters.energy_cap * 100) : INF;
  
  return subopt_par(parameters.sequence, parameters.start_structure, vienna_params, energy_cap, 0, 0, NULL);
}

double estimate_equilibrium(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  int i, everything_in_equilibrium = 0;
  double equilibrium_time;
  
  equilibrium_time = soft_bound_for_population_proportion(eigensystem, parameters, parameters.end_index, parameters.end_time, 1);
  
#ifdef INSANE_DEBUG
  printf("equilibrium_time: %f\n", equilibrium_time);
#endif
  
  while (equilibrium_time >= parameters.start_time && index_in_equilibrium_within_window(eigensystem, parameters, parameters.end_index, equilibrium_time)) {
    equilibrium_time -= parameters.step_size;
  }
  
#ifdef INSANE_DEBUG
  printf("equilibrium_time after backtrack: %f\n", equilibrium_time);
#endif
  
  if (equilibrium_time < parameters.start_time) {
    return NEG_INF(parameters);
  } else {
    while (equilibrium_time <= parameters.end_time) {
      if (parameters.all_subpop_for_eq) {
        everything_in_equilibrium = 1;
        
        for (i = 0; i < eigensystem.length && everything_in_equilibrium; ++i) {
          everything_in_equilibrium = index_in_equilibrium_within_window(eigensystem, parameters, i, equilibrium_time);
        }
      } else {
        everything_in_equilibrium = index_in_equilibrium_within_window(eigensystem, parameters, parameters.end_index, equilibrium_time);
      }
      
      if (everything_in_equilibrium) {
        return equilibrium_time;
      } else {
        equilibrium_time += parameters.step_size;
      }
    }
    
    return POS_INF(parameters);
  }
}

double soft_bound_for_population_proportion(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters, int eigensystem_index, double final_time, int sign) {
  double current_probability, final_probability, current_time, temp_time, last_time = final_time;
  
  if (parameters.soft_bounds) {
    last_time           = final_time;
    current_time        = parameters.start_time + (parameters.end_time - parameters.start_time) / 2.0;
    current_probability = probability_at_logtime(eigensystem, parameters, current_time, parameters.start_index, eigensystem_index);
    final_probability   = probability_at_logtime(eigensystem, parameters, final_time, parameters.start_index, eigensystem_index);
    
    while (fabs(current_time - last_time) > parameters.step_size) {
      current_probability = probability_at_logtime(eigensystem, parameters, current_time, parameters.start_index, eigensystem_index);
      
#ifdef INSANE_DEBUG
      printf("current_time: %+f, last_time: %+f, p(current): %f, p(target): %f\n", current_time, last_time, current_probability, final_probability);
#endif
      
      temp_time = current_time;
      
      if (fabs(current_probability - final_probability) > parameters.delta) {
        current_time += sign * fabs(last_time - current_time) / 2.0;
      } else {
        current_time -= sign * fabs(last_time - current_time) / 2.0;
      }
      
      last_time = temp_time;
    }
    
#ifdef INSANE_DEBUG
    printf("Bounded at %+f within %f of requested delta, %f\n", current_time, fabs(current_probability - final_probability), parameters.delta);
#endif
    
    return round(current_time / parameters.step_size) * parameters.step_size;
  } else {
    return final_time;
  }
}

int index_in_equilibrium_within_window(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters, int eigensystem_index, double logtime) {
  int i;
  double current_proportion, future_proportion;
  
  current_proportion = probability_at_logtime(
                         eigensystem,
                         parameters,
                         logtime,
                         parameters.start_index,
                         eigensystem_index
                       );
                       
  for (i = 1; i < parameters.window_size; ++i) {
    future_proportion = probability_at_logtime(
                          eigensystem,
                          parameters,
                          logtime + parameters.step_size * i,
                          parameters.start_index,
                          eigensystem_index
                        );
                        
#ifdef INSANE_DEBUG
    printf(
      "current_time: %f, window_time: %f, index: %d, abs(p(window) - p(current)): %e, within epsilon?: %s\n",
      logtime,
      logtime + parameters.step_size * i,
      eigensystem_index,
      fabs(current_proportion - future_proportion),
      fabs(current_proportion - future_proportion) < parameters.epsilon ? "yes" : "no"
    );
#endif
    
    if (fabs(current_proportion - future_proportion) > parameters.epsilon) {
      return 0;
    }
  }
  
  return 1;
}

void print_equilibrium(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  double equilibrium_time;
  
  equilibrium_time = estimate_equilibrium(eigensystem, parameters);
  
  if (equilibrium_time == NEG_INF(parameters)) {
    printf("-Infinity\n");
  } else if (equilibrium_time == POS_INF(parameters)) {
    printf("Infinity\n");
  } else {
    printf("%f\n", equilibrium_time);
  }
}

void print_population_proportion(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  double current_time, soft_left_bound, soft_right_bound;
  
  soft_left_bound  = soft_bound_for_population_proportion(eigensystem, parameters, parameters.start_index, parameters.start_time, -1);
  soft_right_bound = soft_bound_for_population_proportion(eigensystem, parameters, parameters.end_index, parameters.end_time, 1);
  
  for (current_time = soft_left_bound; current_time <= soft_right_bound; current_time += parameters.step_size) {
    printf(
      "%+f\t%+.8f\t%+.8f\n",
      current_time,
      probability_at_logtime(eigensystem, parameters, current_time, parameters.start_index, parameters.end_index),
      probability_at_logtime(eigensystem, parameters, current_time, parameters.start_index, parameters.start_index)
    );
  }
}

void print_array(char* title, double* matrix, int length) {
  int i;
  printf("%s\n", title);
  
  for (i = 0; i < length; ++i) {
    printf("%+.4f\n", matrix[i]);
  }
  
  printf("\n");
}


void print_matrix(char* title, double* matrix, int length) {
  int i, j;
  printf("%s\n", title);
  
  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) {
      printf("%+.2e ", E_ROW_ORDER(matrix, i, j, length));
    }
    
    printf("\n");
  }
  
  printf("\n");
}

void print_eigenvalues(const EIGENSYSTEM eigensystem) {
  int i;
  
  for (i = 0; i < eigensystem.length; ++i) {
    printf("%+.8f\n", eigensystem.values[i]);
  }
}
