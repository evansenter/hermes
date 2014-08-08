#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "klp_matrix_initializers.h"
#include "klp_matrix_functions.h"
#include "shared/constants.h"

TRANSITION_MATRIX transpose_matrix(TRANSITION_MATRIX matrix) {
  int i, j;
  TRANSITION_MATRIX transposed_matrix = init_transition_matrix(matrix.row_length, matrix.type);
  
  for (i = 0; i < matrix.row_length; ++i) {
    for (j = 0; j < matrix.row_length; ++j) {
      T_COL_ORDER(transposed_matrix, i, j) = T_ROW_ORDER(matrix, i, j);
    }
  }
  
  free_transition_matrix(matrix);
  return transposed_matrix;
}

TRANSITION_MATRIX inverse(TRANSITION_MATRIX transition_matrix) {
  int i, j, signum;
  gsl_matrix* matrix_to_invert = gsl_matrix_alloc(transition_matrix.row_length, transition_matrix.row_length);
  gsl_matrix* inversion_matrix = gsl_matrix_alloc(transition_matrix.row_length, transition_matrix.row_length);
  gsl_permutation* permutation = gsl_permutation_alloc(transition_matrix.row_length);
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      gsl_matrix_set(matrix_to_invert, i, j, T_ROW_ORDER(transition_matrix, i, j));
    }
  }
  
  gsl_linalg_LU_decomp(matrix_to_invert, permutation, &signum);
  gsl_linalg_LU_invert(matrix_to_invert, permutation, inversion_matrix);
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      T_ROW_ORDER(transition_matrix, i, j) = gsl_matrix_get(inversion_matrix, i, j);
    }
  }
  
  gsl_matrix_free(matrix_to_invert);
  gsl_matrix_free(inversion_matrix);
  gsl_permutation_free(permutation);
  
  return transition_matrix;
}

TRANSITION_MATRIX convert_klp_matrix_to_transition_matrix(KLP_MATRIX* klp_matrix, KLP_PARAMS* klp_params) {
  int resolved;
  double* number_of_adjacent_moves;
  transition_probability probability_function = NULL;
  
  resolved = find_start_and_end_positions_in_klp_matrix(klp_matrix, klp_params);
  
  if (klp_params->max_dist) {
    if (!klp_params->bp_dist) {
      set_bp_dist_from_start_and_end_positions(*klp_matrix, klp_params, resolved);
    }
    
    extend_klp_matrix_to_all_possible_positions(klp_matrix, *klp_params);
    populate_remaining_probabilities_in_klp_matrix(klp_matrix, *klp_params);
    
    if (resolved != 2) {
      find_start_and_end_positions_in_klp_matrix(klp_matrix, klp_params);
    }
  }
  
  number_of_adjacent_moves = populate_number_of_adjacent_moves(*klp_matrix, *klp_params);
  
#ifdef DEBUG
  int i;
  printf("\nFull dataset:\n");
  
  for (i = 0; i < klp_matrix->length; ++i) {
    printf("%d\t%d\t%f\t%d possible move(s)\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i], (int)number_of_adjacent_moves[i]);
  }
  
  printf("\n");
#endif
  
  switch (10 * klp_params->hastings + klp_params->energy_based) {
    case 0:
      probability_function = &transition_rate_from_probabilities;
#ifdef DEBUG
      printf("probability_function: transition_rate_from_probabilities\n");
#endif
      break;
      
    case 1:
      probability_function = &transition_rate_from_energies;
#ifdef DEBUG
      printf("probability_function: transition_rate_from_energies\n");
#endif
      break;
      
    case 10:
      probability_function = &transition_rate_from_probabilities_with_hastings;
#ifdef DEBUG
      printf("probability_function: transition_rate_from_probabilities_with_hastings\n");
#endif
      break;
      
    case 11:
      probability_function = &transition_rate_from_energies_with_hastings;
#ifdef DEBUG
      printf("probability_function: transition_rate_from_energies_with_hastings\n");
#endif
      break;
  }
  
  return populate_transition_matrix_from_stationary_matrix(*klp_matrix, *klp_params, number_of_adjacent_moves, probability_function);
}

int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX* klp_matrix, KLP_PARAMS* klp_params) {
  int i, resolved = 0;
  
  if (klp_params->start_state == -1) {
    for (i = 0; i < klp_matrix->length && klp_params->start_state == -1; ++i) {
      if (klp_matrix->k[i] == 0) {
        klp_params->start_state = i;
        resolved++;
      }
    }
  } else {
    resolved++;
  }
  
  if (klp_params->end_state == -1) {
    for (i = 0; i < klp_matrix->length && klp_params->end_state == -1; ++i) {
      if (klp_matrix->l[i] == 0) {
        klp_params->end_state = i;
        resolved++;
      }
    }
  } else {
    resolved++;
  }
  
#ifdef DEBUG
  printf("\nstart_index:\t%d\n", klp_params->start_state);
  printf("end_index:\t%d\n", klp_params->end_state);
  printf("resolved:\t%d\n", resolved);
#endif
  return resolved;
}

void set_bp_dist_from_start_and_end_positions(const KLP_MATRIX klp_matrix, KLP_PARAMS* klp_params, int resolved) {
  int distance_from_start, distance_from_end;
  
  distance_from_start = distance_from_end = -1;
  
  if (klp_params->start_state > 0) {
    distance_from_start = klp_matrix.l[klp_params->start_state];
  }
  
  if (klp_params->end_state > 0) {
    distance_from_end = klp_matrix.k[klp_params->end_state];
  }
  
  if (distance_from_start == distance_from_end && resolved) {
    klp_params->bp_dist = distance_from_start;
  } else if (distance_from_start >= 0 && distance_from_end == -1) {
    klp_params->bp_dist = distance_from_start;
  } else if (distance_from_end >= 0 && distance_from_start == -1) {
    klp_params->bp_dist = distance_from_end;
  } else {
    fprintf(stderr, "Can't infer the input structure distances for the energy grid. We found (0, %d) and (%d, 0). Consider using the -d flag to manually set the base pair distance between the two structures.\n", distance_from_end, distance_from_start);
    printf("-3\n");
    exit(0);
  }
  
#ifdef DEBUG
  printf("bp_dist:\t%d\n", klp_params->bp_dist);
#endif
}

void extend_klp_matrix_to_all_possible_positions(KLP_MATRIX* klp_matrix, const KLP_PARAMS klp_params) {
  int i, j, m, position_in_input_data, pointer, valid_positions = 0;
  
#ifdef DEBUG
  printf("\nAccessible positions (top-left is [0, 0]):\n");
#endif
  
  for (i = 0; i <= klp_params.max_dist; ++i) {
    for (j = 0; j <= klp_params.max_dist; ++j) {
      if (
        i + j >= klp_params.bp_dist &&
        i + klp_params.bp_dist >= j &&
        j + klp_params.bp_dist >= i &&
        (i + j) % 2 == klp_params.bp_dist % 2
      ) {
#ifdef DEBUG
        position_in_input_data = -1;
        
        for (m = 0; m < klp_matrix->length && position_in_input_data == -1; ++m) {
          if (klp_matrix->k[m] == i && klp_matrix->l[m] == j) {
            position_in_input_data = m;
          }
        }
        
        printf(position_in_input_data == -1 ? "X" : "O");
#endif
        valid_positions++;
      } else {
#ifdef DEBUG
        printf(" ");
#endif
      }
    }
    
#ifdef DEBUG
    printf("\n");
#endif
  }
  
  klp_matrix->k = realloc(klp_matrix->k, valid_positions * sizeof(int));
  klp_matrix->l = realloc(klp_matrix->l, valid_positions * sizeof(int));
  klp_matrix->p = realloc(klp_matrix->p, valid_positions * sizeof(double));
  
  pointer = klp_matrix->length;
  
#ifdef DEBUG
  printf("\nInput dataset:\n");
  
  for (i = 0; i < klp_matrix->length; ++i) {
    printf("%d\t%d\t%d\t%f\n", i, klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i]);
  }
  
#endif
  
  for (i = 0; i <= klp_params.max_dist; ++i) {
    for (j = 0; j <= klp_params.max_dist; ++j) {
      if (
        i + j >= klp_params.bp_dist &&
        i + klp_params.bp_dist >= j &&
        j + klp_params.bp_dist >= i &&
        (i + j) % 2 == klp_params.bp_dist % 2
      ) {
        position_in_input_data = -1;
        
        for (m = 0; m < klp_matrix->length && position_in_input_data == -1; ++m) {
          if (klp_matrix->k[m] == i && klp_matrix->l[m] == j) {
            position_in_input_data = m;
          }
        }
        
        if (position_in_input_data < 0) {
          klp_matrix->k[pointer] = i;
          klp_matrix->l[pointer] = j;
          klp_matrix->p[pointer] = 0.;
          pointer++;
        }
      }
    }
  }
  
  klp_matrix->length = valid_positions;
}

void populate_remaining_probabilities_in_klp_matrix(KLP_MATRIX* klp_matrix, const KLP_PARAMS klp_params) {
  int i;
  double epsilon_per_cell;
  
  if (klp_params.epsilon) {
    // Extend the energy grid by adding an epsilon value to all 0-probability positions.
    epsilon_per_cell = klp_params.epsilon / klp_matrix->length;
    
    for (i = 0; i < klp_matrix->length; ++i) {
      if (klp_matrix->p[i] > 0) {
        klp_matrix->p[i] = (klp_matrix->p[i] + epsilon_per_cell) / (1. + klp_params.epsilon);
      } else {
        klp_matrix->p[i] = epsilon_per_cell / (1. + klp_params.epsilon);
      }
    }
  }
}

double* populate_number_of_adjacent_moves(const KLP_MATRIX klp_matrix, const KLP_PARAMS klp_params) {
  int i;
  double* number_of_adjacent_moves;
  
  number_of_adjacent_moves = malloc(klp_matrix.length * sizeof(double));
  
  for (i = 0; i < klp_matrix.length; ++i) {
    number_of_adjacent_moves[i] = RUN_TYPE(klp_params.run_type, DIAG_MOVES_ONLY_FLAG) ? (double)number_of_permissible_single_bp_moves(klp_matrix, i) : (double)(klp_matrix.length - 1);
  }
  
  return number_of_adjacent_moves;
}

int number_of_permissible_single_bp_moves(const KLP_MATRIX klp_matrix, int i) {
  int j, x, y, a, b, num_moves = 0;
  
  x = klp_matrix.k[i];
  y = klp_matrix.l[i];
  
  for (j = 0; j < klp_matrix.length; ++j) {
    a = klp_matrix.k[j];
    b = klp_matrix.l[j];
    
    if (
      // Because N(x, y) is restricted to entries in *k and *l, we *assume* the input data satisfies the triangle inequality and bounds.
      (int)abs(x - a) == 1 && (int)abs(y - b) == 1
    ) {
      num_moves++;
    }
  }
  
  return num_moves;
}

TRANSITION_MATRIX populate_transition_matrix_from_stationary_matrix(const KLP_MATRIX klp_matrix, const KLP_PARAMS klp_params, const double* number_of_adjacent_moves, transition_probability probability_function) {
  int i, j;
  double row_sum;
  TRANSITION_MATRIX transition_matrix;
  
  transition_matrix = init_transition_matrix(klp_matrix.length, MATRIX_TYPE(klp_params.rate_matrix));
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    row_sum = 0.;
    
    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != j) {
        if (RUN_TYPE(klp_params.run_type, FULLY_CONNECTED_FLAG) || (RUN_TYPE(klp_params.run_type, DIAG_MOVES_ONLY_FLAG) && ONE_BP_MOVE(i, j))) {
          if (NONZERO_TO_NONZERO_PROB(i, j)) {
            T_ROW_ORDER(transition_matrix, i, j) = \
                                                   probability_function(klp_matrix, number_of_adjacent_moves, i, j, klp_params.rate_matrix);
          }
        }
        
        row_sum += T_ROW_ORDER(transition_matrix, i, j);
      }
    }
    
    T_ROW_ORDER(transition_matrix, i, i) = klp_params.rate_matrix ? -row_sum : 1 - row_sum;
  }
  
  return transition_matrix;
}

double transition_rate_from_probabilities(const KLP_MATRIX klp_matrix, const double* number_of_adjacent_moves, int i, int j, short rate_matrix) {
  if (rate_matrix) {
    return MIN(1., klp_matrix.p[j] / klp_matrix.p[i]);
  } else {
    return MIN(1., klp_matrix.p[j] / klp_matrix.p[i]) / number_of_adjacent_moves[i];
  }
}

double transition_rate_from_energies(const KLP_MATRIX klp_matrix, const double* number_of_adjacent_moves, int i, int j, short rate_matrix) {
  if (rate_matrix) {
    return MIN(1., exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT));
  } else {
    return MIN(1., exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT)) / number_of_adjacent_moves[i];
  }
}

double transition_rate_from_probabilities_with_hastings(const KLP_MATRIX klp_matrix, const double* number_of_adjacent_moves, int i, int j, short rate_matrix) {
  if (rate_matrix) {
    return MIN(1., (number_of_adjacent_moves[i] / number_of_adjacent_moves[j]) * (klp_matrix.p[j] / klp_matrix.p[i]));
  } else {
    return MIN(1., (number_of_adjacent_moves[i] / number_of_adjacent_moves[j]) * (klp_matrix.p[j] / klp_matrix.p[i])) / number_of_adjacent_moves[i];
  }
}

double transition_rate_from_energies_with_hastings(const KLP_MATRIX klp_matrix, const double* number_of_adjacent_moves, int i, int j, short rate_matrix) {
  if (rate_matrix) {
    return MIN(1., (number_of_adjacent_moves[i] / number_of_adjacent_moves[j]) * exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT));
  } else {
    return MIN(1., (number_of_adjacent_moves[i] / number_of_adjacent_moves[j]) * exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT)) / number_of_adjacent_moves[i];
  }
}
