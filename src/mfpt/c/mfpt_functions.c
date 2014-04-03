#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "constants.h"
#include "initializers.h"
#include "functions.h"

#define ONE_BP_MOVE(i, j) ((int)abs(klp_matrix.k[(i)] - klp_matrix.k[(j)]) == 1 && (int)abs(klp_matrix.l[(i)] - klp_matrix.l[(j)]) == 1)
#define NONZERO_TO_NONZERO_PROB(i, j) (klp_matrix.p[(i)] > 0 && klp_matrix.p[(j)] > 0)

TRANSITION_MATRIX convert_klp_matrix_to_transition_matrix(KLP_MATRIX* klp_matrix, MFPT_PARAMS* parameters) {
  int resolved;
  double* number_of_adjacent_moves;
  transition_probability probability_function = NULL;

  resolved = find_start_and_end_positions_in_klp_matrix(klp_matrix, parameters);

  if (parameters->max_dist) {
    if (!parameters->bp_dist) {
      set_bp_dist_from_start_and_end_positions(*klp_matrix, parameters, resolved);
    }

    extend_klp_matrix_to_all_possible_positions(klp_matrix, *parameters);
    populate_remaining_probabilities_in_klp_matrix(klp_matrix, *parameters);

    if (resolved != 2) {
      find_start_and_end_positions_in_klp_matrix(klp_matrix, parameters);
    }
  }

  number_of_adjacent_moves = populate_number_of_adjacent_moves(*klp_matrix, *parameters);

#ifdef DEBUG
  int i;
  printf("\nFull dataset:\n");

  for (i = 0; i < klp_matrix->length; ++i) {
    printf("%d\t%d\t%f\t%d possible move(s)\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i], (int)number_of_adjacent_moves[i]);
  }

  printf("\n");
#endif

  switch (10 * parameters->hastings + parameters->energy_based) {
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

  return populate_transition_matrix_from_stationary_matrix(*klp_matrix, *parameters, number_of_adjacent_moves, probability_function);
}

double compute_mfpt(const TRANSITION_MATRIX transition_matrix, const MFPT_PARAMS parameters) {
  int i, j, x, y, start_pointer;
  double mfpt_from_start;
  TRANSITION_MATRIX inversion_matrix;
  double* mfpt;

  if (parameters.start_state < 0 || parameters.end_state < 0) {
    if (parameters.start_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
#endif
      return -1;
    }

    if (parameters.end_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the stopping state.\n");
#endif
      return -2;
    }
  }

  // If start_index > end_index, we need to shift to the left by one because the end_index row / column is being removed.
  start_pointer    = parameters.start_state - (parameters.start_state > parameters.end_state ? 1 : 0);
  inversion_matrix = init_transition_matrix(transition_matrix.row_length - 1, transition_matrix.type);
  mfpt             = calloc(inversion_matrix.row_length, sizeof(double));

  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != parameters.end_state && j != parameters.end_state) {
        x = (i > parameters.end_state ? i - 1 : i);
        y = (j > parameters.end_state ? j - 1 : j);
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (end_index) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversion_matrix and i, j are used for indexing into transition_matrix.
        ROW_ORDER(inversion_matrix, x, y) = \
            (i == j ? 1 - ROW_ORDER(transition_matrix, i, j) : -ROW_ORDER(transition_matrix, i, j));
      }
    }
  }

  inversion_matrix = parameters.pseudoinverse ? pseudoinverse(inversion_matrix) : inverse(inversion_matrix);

  for (i = 0; i < inversion_matrix.row_length; ++i) {
    for (j = 0; j < inversion_matrix.row_length; ++j) {
      mfpt[i] += ROW_ORDER(inversion_matrix, i, j);
    }

    if (parameters.all_mfpt) {
      // The business with this i < end_index stuff is inorder to ensure that the output MFPT indices are representative of the input data, since we reduce the dimension of the matrix by 1.
      printf("%d\t%+.8f\n", i < parameters.end_state ? i : i + 1, mfpt[i]);
    }
  }

  mfpt_from_start = mfpt[start_pointer];

  free(mfpt);
  free_transition_matrix(inversion_matrix);
  return mfpt_from_start;
}

TRANSITION_MATRIX inverse(TRANSITION_MATRIX matrix_to_invert) {
  // http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
  double* a = matrix_to_invert.matrix;
  int size  = matrix_to_invert.row_length;

  // -----------------------------------------------------
  // Begin LAPACK calls
  // -----------------------------------------------------
  int* ipiv    = malloc((size + 1) * sizeof(int));
  int lwork    = size * size;
  double* work = malloc(lwork * sizeof(double));
  int info;
  dgetrf_(&size, &size, a, &size, ipiv, &info);
#ifdef DEBUG
  printf("dgetrf_(&size, &size, a, &size, ipiv, &info)\ninfo:\t%d\n", info);
#endif
  dgetri_(&size, a, &size, ipiv, work, &lwork, &info);
#ifdef DEBUG
  printf("dgetri_(&size, a, &size, ipiv, work, &lwork, &info)\ninfo:\t%d\n", info);
#endif
  free(ipiv);
  free(work);
  // -----------------------------------------------------
  // End LAPACK calls
  // -----------------------------------------------------

  return matrix_to_invert;
}

TRANSITION_MATRIX pseudoinverse(TRANSITION_MATRIX matrix_to_invert) {
  // Least-squares fit solution to B - Ax, where (in this case) A is square and B is the identity matrix.
  double* a = matrix_to_invert.matrix;
  int size  = matrix_to_invert.row_length;

  // -----------------------------------------------------
  // Begin LAPACK calls
  // -----------------------------------------------------
  char trans;
  int i, m, n, nrhs, lda, ldb, lwork, info;
  trans = 'N';
  m     = size;
  n     = m;
  nrhs  = m;
  lda   = m;
  ldb   = m;
  double* b = calloc(ldb * nrhs, sizeof(double));

  for (i = 0; i < ldb; ++i) {
    b[i * nrhs + i] = 1.;
  }

#ifdef SUPER_HEAVY_DEBUG
  printf("dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info)\n\n");
#endif
  lwork        = -1;
  double* work = malloc(MAX(1, lwork) * sizeof(double));
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  lwork = (int)work[0];
#ifdef SUPER_HEAVY_DEBUG
  printf("workspace query (lwork):\t%d\n", lwork);
#endif
  free(work);
  work = malloc(MAX(1, lwork) * sizeof(double));
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
#ifdef DEBUG
  printf("info:\t%d\n", info);
#endif
  free(work);
  // -----------------------------------------------------
  // End LAPACK calls.
  // -----------------------------------------------------

  TRANSITION_MATRIX inverted_matrix = {
    .matrix     = b,
    .row_length = matrix_to_invert.row_length,
    .type       = matrix_to_invert.type
  };

  free_transition_matrix(matrix_to_invert);

  return inverted_matrix;
}

int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX* klp_matrix, MFPT_PARAMS* parameters) {
  int i, resolved = 0;

  if (parameters->start_state == -1) {
    for (i = 0; i < klp_matrix->length && parameters->start_state == -1; ++i) {
      if (klp_matrix->k[i] == 0) {
        parameters->start_state = i;
        resolved++;
      }
    }
  } else {
    resolved++;
  }

  if (parameters->end_state == -1) {
    for (i = 0; i < klp_matrix->length && parameters->end_state == -1; ++i) {
      if (klp_matrix->l[i] == 0) {
        parameters->end_state = i;
        resolved++;
      }
    }
  } else {
    resolved++;
  }

#ifdef DEBUG
  printf("\nstart_index:\t%d\n", parameters->start_state);
  printf("end_index:\t%d\n", parameters->end_state);
  printf("resolved:\t%d\n", resolved);
#endif
  return resolved;
}

void set_bp_dist_from_start_and_end_positions(const KLP_MATRIX klp_matrix, MFPT_PARAMS* parameters, int resolved) {
  int distance_from_start, distance_from_end;

  distance_from_start = distance_from_end = -1;

  if (parameters->start_state > 0) {
    distance_from_start = klp_matrix.l[parameters->start_state];
  }

  if (parameters->end_state > 0) {
    distance_from_end = klp_matrix.k[parameters->end_state];
  }

  if (distance_from_start == distance_from_end && resolved) {
    parameters->bp_dist = distance_from_start;
  } else if (distance_from_start >= 0 && distance_from_end == -1) {
    parameters->bp_dist = distance_from_start;
  } else if (distance_from_end >= 0 && distance_from_start == -1) {
    parameters->bp_dist = distance_from_end;
  } else {
    fprintf(stderr, "Can't infer the input structure distances for the energy grid. We found (0, %d) and (%d, 0). Consider using the -d flag to manually set the base pair distance between the two structures.\n", distance_from_end, distance_from_start);
    printf("-3\n");
    exit(0);
  }

#ifdef DEBUG
  printf("bp_dist:\t%d\n", parameters->bp_dist);
#endif
}

void extend_klp_matrix_to_all_possible_positions(KLP_MATRIX* klp_matrix, const MFPT_PARAMS parameters) {
  int i, j, m, position_in_input_data, pointer, valid_positions = 0;

#ifdef DEBUG
  printf("\nAccessible positions (top-left is [0, 0]):\n");
#endif

  for (i = 0; i <= parameters.max_dist; ++i) {
    for (j = 0; j <= parameters.max_dist; ++j) {
      if (
        i + j >= parameters.bp_dist &&
        i + parameters.bp_dist >= j &&
        j + parameters.bp_dist >= i &&
        (i + j) % 2 == parameters.bp_dist % 2
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

  for (i = 0; i <= parameters.max_dist; ++i) {
    for (j = 0; j <= parameters.max_dist; ++j) {
      if (
        i + j >= parameters.bp_dist &&
        i + parameters.bp_dist >= j &&
        j + parameters.bp_dist >= i &&
        (i + j) % 2 == parameters.bp_dist % 2
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

void populate_remaining_probabilities_in_klp_matrix(KLP_MATRIX* klp_matrix, const MFPT_PARAMS parameters) {
  int i;
  double epsilon_per_cell;

  if (parameters.epsilon) {
    // Extend the energy grid by adding an epsilon value to all 0-probability positions.
    epsilon_per_cell = parameters.epsilon / klp_matrix->length;

    for (i = 0; i < klp_matrix->length; ++i) {
      if (klp_matrix->p[i] > 0) {
        klp_matrix->p[i] = (klp_matrix->p[i] + epsilon_per_cell) / (1. + parameters.epsilon);
      } else {
        klp_matrix->p[i] = epsilon_per_cell / (1. + parameters.epsilon);
      }
    }
  }
}

double* populate_number_of_adjacent_moves(const KLP_MATRIX klp_matrix, const MFPT_PARAMS parameters) {
  int i;
  double* number_of_adjacent_moves;

  number_of_adjacent_moves = malloc(klp_matrix.length * sizeof(double));

  for (i = 0; i < klp_matrix.length; ++i) {
    number_of_adjacent_moves[i] = RUN_TYPE(parameters, DIAG_MOVES_ONLY_FLAG) ? (double)number_of_permissible_single_bp_moves(klp_matrix, i) : (double)(klp_matrix.length - 1);
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

TRANSITION_MATRIX populate_transition_matrix_from_stationary_matrix(const KLP_MATRIX klp_matrix, const MFPT_PARAMS parameters, const double* number_of_adjacent_moves, transition_probability probability_function) {
  int i, j;
  double row_sum;
  TRANSITION_MATRIX transition_matrix;

  transition_matrix = init_transition_matrix(klp_matrix.length, MATRIX_TYPE(parameters));

  for (i = 0; i < transition_matrix.row_length; ++i) {
    row_sum = 0.;

    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != j) {
        if (RUN_TYPE(parameters, FULLY_CONNECTED_FLAG) || (RUN_TYPE(parameters, DIAG_MOVES_ONLY_FLAG) && ONE_BP_MOVE(i, j))) {
          if (NONZERO_TO_NONZERO_PROB(i, j)) {
            ROW_ORDER(transition_matrix, i, j) = \
                probability_function(klp_matrix, number_of_adjacent_moves, i, j, parameters.rate_matrix);
          }
        }

        row_sum += ROW_ORDER(transition_matrix, i, j);
      }
    }

    ROW_ORDER(transition_matrix, i, i) = parameters.rate_matrix ? -row_sum : 1 - row_sum;
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
