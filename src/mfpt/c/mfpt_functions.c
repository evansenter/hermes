#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "mfpt_functions.h"
#include "shared/libklp_matrix_header.h"
#include "shared/constants.h"

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

double compute_mfpt(const TRANSITION_MATRIX transition_matrix, const KLP_PARAMS klp_params, const MFPT_PARAMS parameters) {
  int i, j, x, y, start_pointer;
  double mfpt_from_start;
  TRANSITION_MATRIX inversion_matrix;
  double* mfpt;
  
  if (klp_params.start_state < 0 || klp_params.end_state < 0) {
    if (klp_params.start_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
#endif
      return -1;
    }
    
    if (klp_params.end_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the stopping state.\n");
#endif
      return -2;
    }
  }
  
  // If start_index > end_index, we need to shift to the left by one because the end_index row / column is being removed.
  start_pointer    = klp_params.start_state - (klp_params.start_state > klp_params.end_state ? 1 : 0);
  inversion_matrix = init_transition_matrix(transition_matrix.row_length - 1, transition_matrix.type);
  mfpt             = calloc(inversion_matrix.row_length, sizeof(double));
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != klp_params.end_state && j != klp_params.end_state) {
        x = (i > klp_params.end_state ? i - 1 : i);
        y = (j > klp_params.end_state ? j - 1 : j);
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (end_index) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversion_matrix and i, j are used for indexing into transition_matrix.
        T_ROW_ORDER(inversion_matrix, x, y) = \
                                              (i == j ? 1 - T_ROW_ORDER(transition_matrix, i, j) : -T_ROW_ORDER(transition_matrix, i, j));
      }
    }
  }
  
  inversion_matrix = inverse(inversion_matrix);
  
  for (i = 0; i < inversion_matrix.row_length; ++i) {
    for (j = 0; j < inversion_matrix.row_length; ++j) {
      mfpt[i] += T_ROW_ORDER(inversion_matrix, i, j);
    }
    
    if (parameters.all_mfpt) {
      // The business with this i < end_index stuff is inorder to ensure that the output MFPT indices are representative of the input data, since we reduce the dimension of the matrix by 1.
      printf("%d\t%+.8f\n", i < klp_params.end_state ? i : i + 1, mfpt[i]);
    }
  }
  
  mfpt_from_start = mfpt[start_pointer];
  
  free(mfpt);
  free_transition_matrix(inversion_matrix);
  return mfpt_from_start;
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
    number_of_adjacent_moves[i] = RUN_TYPE(klp_params, DIAG_MOVES_ONLY_FLAG) ? (double)number_of_permissible_single_bp_moves(klp_matrix, i) : (double)(klp_matrix.length - 1);
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
  
  transition_matrix = init_transition_matrix(klp_matrix.length, MATRIX_TYPE(klp_params));
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    row_sum = 0.;
    
    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != j) {
        if (RUN_TYPE(klp_params, FULLY_CONNECTED_FLAG) || (RUN_TYPE(klp_params, DIAG_MOVES_ONLY_FLAG) && ONE_BP_MOVE(i, j))) {
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
