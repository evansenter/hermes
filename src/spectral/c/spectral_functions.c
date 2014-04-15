#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "shared/libmfpt_header.h"
#include "shared/libtpl_header.h"
#include "constants.h"
#include "initializers.h"
#include "functions.h"

EIGENSYSTEM eigensystem_from_row_ordered_transition_matrix(TRANSITION_MATRIX row_transition_matrix) {
  TRANSITION_MATRIX column_transition_matrix;
  EIGENSYSTEM eigensystem;
  
  column_transition_matrix = transpose_matrix(row_transition_matrix);
  eigensystem              = convert_transition_matrix_to_eigenvectors(column_transition_matrix);
  invert_matrix(&eigensystem);
  
  return eigensystem;
}

void population_proportion_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
  EIGENSYSTEM eigensystem;
  
  eigensystem = eigensystem_from_row_ordered_transition_matrix(row_transition_matrix);
  print_population_proportion(eigensystem, parameters);
}

void equilibrium_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
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

double probability_at_time(const EIGENSYSTEM eigensystem, double timepoint, int start_index, int target_index) {
  // This function is hard-wired to only consider the kinetics for folding from a distribution where p_{0}(start_index) == 1.
  int i;
  double cumulative_probability = 0;
  
  for (i = 0; i < eigensystem.length; ++i) {
    cumulative_probability +=
      E_ROW_ORDER(eigensystem.vectors, target_index, i, eigensystem.length) *
      E_ROW_ORDER(eigensystem.inverse_vectors, i, start_index, eigensystem.length) *
      exp(eigensystem.values[i] * timepoint);
  }
  
  return cumulative_probability;
}

void find_key_structure_indices_in_structure_list(SPECTRAL_PARAMS* parameters, const SOLUTION* all_structures, int num_structures, char* empty_str, char* mfe_str) {
  int i;
  
  for (i = 0; i < num_structures; ++i) {
    if (parameters->start_index == -1) {
      if (parameters->start_structure != NULL) {
        if (!strcmp(parameters->start_structure, all_structures[i].structure)) {
          parameters->start_index = i;
        }
      } else {
        if (!strcmp(empty_str, all_structures[i].structure)) {
          parameters->start_structure = all_structures[i].structure;
          parameters->start_index     = i;
        }
      }
    }
    
    if (parameters->end_index == -1) {
      if (parameters->end_structure != NULL) {
        if (!strcmp(parameters->end_structure, all_structures[i].structure)) {
          parameters->end_index = i;
        }
      } else {
        if (!strcmp(mfe_str, all_structures[i].structure)) {
          parameters->end_structure = all_structures[i].structure;
          parameters->end_index     = i;
        }
      }
    }
  }
}

void serialize_eigensystem(const EIGENSYSTEM eigensystem, const SPECTRAL_PARAMS parameters) {
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

EIGENSYSTEM deserialize_eigensystem(const SPECTRAL_PARAMS parameters) {
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

long double estimate_equilibrium(const EIGENSYSTEM eigensystem, const SPECTRAL_PARAMS parameters) {
  int i, j, starting_index, num_points, equilibrium_count = 0;
  double current_proportion, future_proportion, epsilon = parameters.equilibrium;
  
  num_points     = (int)((parameters.end_time - parameters.start_time) / parameters.step_size);
  starting_index = estimate_starting_index_to_scan_for_equilibrium(num_points, eigensystem, parameters);
  
  for (i = starting_index; i < num_points - WINDOW_SIZE + 1 && equilibrium_count < eigensystem.length; ++i) {
    equilibrium_count = 0;
    
    for (j = 0; j < eigensystem.length; ++j) {
      current_proportion = probability_at_time(
                             eigensystem,
                             pow(10, parameters.start_time + parameters.step_size * i),
                             parameters.start_index,
                             j
                           );
                           
      future_proportion = probability_at_time(
                            eigensystem,
                            pow(10, parameters.start_time + parameters.step_size * (i + WINDOW_SIZE - 1)),
                            parameters.start_index,
                            j
                          );
                          
#ifdef INSANE_DEBUG
      printf("%d\t%d\t%f\t%f\t%e\n", i, j, current_proportion, future_proportion, fabs(current_proportion - future_proportion));
#endif
      
      if (fabs(current_proportion - future_proportion) < epsilon) {
        equilibrium_count++;
      }
    }
  }
  
#ifdef DEBUG
  printf("eq_index:\t%d\n", i);
#endif
  
  if (i == num_points - WINDOW_SIZE + 1) {
    return -1;
  } else {
    return pow(10, parameters.start_time + parameters.step_size * i);
  }
}

int estimate_starting_index_to_scan_for_equilibrium(int num_points, const EIGENSYSTEM eigensystem, const SPECTRAL_PARAMS parameters) {
  int i = 0, starting_index = 0;
  double str_1_start, str_1_current, str_1_end, str_2_start, str_2_current, str_2_end, epsilon = parameters.equilibrium;
  
  str_1_start = probability_at_time(eigensystem, pow(10, parameters.start_time), parameters.start_index, parameters.start_index);
  str_1_end   = probability_at_time(eigensystem, pow(10, parameters.end_time),   parameters.start_index, parameters.start_index);
  
  str_2_start = probability_at_time(eigensystem, pow(10, parameters.start_time), parameters.start_index, parameters.end_index);
  str_2_end   = probability_at_time(eigensystem, pow(10, parameters.end_time),   parameters.start_index, parameters.end_index);
  
#ifdef DEBUG
  printf("epsilon:\t%.2e\n", epsilon);
  printf("str_1_start:\t%f\n", str_1_start);
  printf("str_1_end:\t%f\n", str_1_end);
  printf("str_2_start:\t%f\n", str_2_start);
  printf("str_2_end:\t%f\n", str_2_end);
#endif
  
  if (fabs(str_1_start - str_1_end) > epsilon || fabs(str_2_start - str_2_end) > epsilon) {
    do {
      str_1_current = probability_at_time(
                        eigensystem,
                        pow(10, parameters.start_time + parameters.step_size * i),
                        parameters.start_index,
                        parameters.start_index
                      );
                      
      str_2_current = probability_at_time(
                        eigensystem,
                        pow(10, parameters.start_time + parameters.step_size * i),
                        parameters.start_index,
                        parameters.end_index
                      );
                      
#ifdef INSANE_DEBUG
      printf("%d\t%f\t%f\n", i, fabs(str_1_current - str_1_end), fabs(str_2_current - str_2_end));
#endif
      
      starting_index = i++;
    } while (i < num_points - WINDOW_SIZE && (fabs(str_1_current - str_1_end) > epsilon || fabs(str_2_current - str_2_end) > epsilon));
  }
  
#ifdef DEBUG
  printf("starting_index:\t%d\n", starting_index);
  printf("num_points:\t%d\n", num_points);
#endif
  
  return starting_index;
}

void print_equilibrium(const EIGENSYSTEM eigensystem, const SPECTRAL_PARAMS parameters) {
  long double equilibrium_time;
  
  equilibrium_time = estimate_equilibrium(eigensystem, parameters);
  
  if (equilibrium_time > 0) {
    printf("%f\n", log(equilibrium_time) / log(10.));
  } else {
    printf("Infinity\n");
  }
}

void print_population_proportion(const EIGENSYSTEM eigensystem, const SPECTRAL_PARAMS parameters) {
  double step_counter;
  
  for (step_counter = parameters.start_time; step_counter <= parameters.end_time; step_counter += parameters.step_size) {
    printf(
      "%+f\t%+.8f\t%+.8f\n",
      step_counter,
      probability_at_time(eigensystem, pow(10, step_counter), parameters.start_index, parameters.end_index),
      probability_at_time(eigensystem, pow(10, step_counter), parameters.start_index, parameters.start_index)
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
