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

void equilibrium_from_row_ordered_transition_matrix(const POPULATION_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
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

double probability_at_time(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters, double timepoint, int start_index, int target_index) {
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
      exp(eigensystem.values[i] * timepoint);
  }

  return cumulative_probability;
}

void find_key_structure_indices_in_structure_list(POPULATION_PARAMS* parameters, const SOLUTION* all_structures, int num_structures, char* empty_str, char* mfe_str) {
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

long double estimate_equilibrium(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  int i, j, starting_index, num_points, everything_in_equilibrium = 0;
  double start_probability, end_probability, epsilon = parameters.equilibrium;

  // estimate_starting_index_to_scan_for_equilibrium does a forward / backward scan to pick the first position to use a sliding window looking for equilibrium.
  // This is necessary because if the timespan is extended too far to the left, the algorithm would erroneously claim equilibrium had been achieved when in fact,
  // the system hadn't had enough time to start making transitions yet. (That corresponds to the left-hand side of the sigmoidal population occupancy curve).
  num_points     = (int)((parameters.end_time - parameters.start_time) / parameters.step_size);
  starting_index = estimate_starting_index_to_scan_for_equilibrium(num_points, eigensystem, parameters);

  // starting_index < 0 means that for the delta provided, we can't find a starting position to scan for equilibrium. This usually means the timespan isn't
  // wide enough.
  if (starting_index < 0) {
    if (parameters.all_subpop_for_eq) {
      everything_in_equilibrium = 1;

      // This just makes sure that if the user is looking at all subpopulations, they all are in equilibrium for the time requested.
      for (i = 0; i < eigensystem.length && everything_in_equilibrium; ++i) {
        start_probability = probability_at_time(eigensystem, parameters, pow(10, parameters.start_time), parameters.start_index, i);
        end_probability   = probability_at_time(eigensystem, parameters, pow(10, parameters.end_time),   parameters.start_index, i);

        everything_in_equilibrium = fabs(start_probability - end_probability) < epsilon;
      }
    } else {
      everything_in_equilibrium = 1;
    }

    // If everything is in equilibrium for the entire time range, we can't compute the equilibrium time and return -Infinity as a result.
    if (everything_in_equilibrium) {
      return NEG_INF;
    }
  }

  // This is a greedy approach. The outer loop assumes that we are not in eq. and walks along time. The inner loop assumes we are in eq. and walks along the
  // subcollections represented in the eigensystem. The first time a failure condition is met in the inner loop (i.e. we're not in equilibrium for that window /
  // subcollection according to the function), the inner loop early terminates after having reflagged the everything_in_equilibrium as false (the returned)
  // value from the function. The only success condition is every evaluation of the inner loop reassigning everything_in_equilibrium to 1, preserving the bit
  // and exiting the outer loop. A goto / (do / while) probably would have worked here as well but I don't think those read as well. Then again, this comment
  // is pretty long so this probably doesn't read that well either.
  for (i = starting_index; i < num_points - parameters.window_size + 1 && !everything_in_equilibrium; ++i) {
    if (parameters.all_subpop_for_eq) {
      everything_in_equilibrium = 1;

      for (j = 0; j < eigensystem.length && everything_in_equilibrium; ++j) {
        everything_in_equilibrium = is_index_in_equilibrium_within_window_position(eigensystem, parameters, j, i);
      }
    } else {
      everything_in_equilibrium = is_index_in_equilibrium_within_window_position(eigensystem, parameters, parameters.end_index, i);
    }
  }

#ifdef DEBUG
  printf("eq_index:\t%d\n", i);
#endif

  if (everything_in_equilibrium) {
    return pow(10, parameters.start_time + parameters.step_size * i);
  } else {
    return POS_INF;
  }
}

int estimate_starting_index_to_scan_for_equilibrium(int num_points, const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  // This scans for a good starting position to estimate equilibrium in a "forward / backward" fashion, described below.
  int i = 0, starting_index = 0;
  double str_1_start, str_1_current, str_1_end, str_2_start, str_2_current, str_2_end, delta = parameters.delta, epsilon = parameters.equilibrium;

  if (!delta) {
    delta = epsilon;
  }

  str_1_start = probability_at_time(eigensystem, parameters, pow(10, parameters.start_time), parameters.start_index, parameters.start_index);
  str_1_end   = probability_at_time(eigensystem, parameters, pow(10, parameters.end_time),   parameters.start_index, parameters.start_index);

  str_2_start = probability_at_time(eigensystem, parameters, pow(10, parameters.start_time), parameters.start_index, parameters.end_index);
  str_2_end   = probability_at_time(eigensystem, parameters, pow(10, parameters.end_time),   parameters.start_index, parameters.end_index);

#ifdef DEBUG
  printf("epsilon:\t%.2e\n",   epsilon);
  printf("str_1_start:\t%f\n", str_1_start);
  printf("str_1_end:\t%f\n",   str_1_end);
  printf("str_2_start:\t%f\n", str_2_start);
  printf("str_2_end:\t%f\n",   str_2_end);
#endif

  // The forward part of selecting a starting index first requires ensuring that a solution exists, namely that for the user-requested timespan, the system
  // (at least for either the population proportion of the starting state or ending state) changes by more than delta overall. If delta is not provided,
  // epsilon is used instead.
  if (fabs(str_1_start - str_1_end) > delta || fabs(str_2_start - str_2_end) > delta) {
    do {
      str_1_current = probability_at_time(
                        eigensystem,
                        parameters,
                        pow(10, parameters.start_time + parameters.step_size * i),
                        parameters.start_index,
                        parameters.start_index
                      );

      str_2_current = probability_at_time(
                        eigensystem,
                        parameters,
                        pow(10, parameters.start_time + parameters.step_size * i),
                        parameters.start_index,
                        parameters.end_index
                      );

#ifdef INSANE_DEBUG
      printf("%d\t%f\t%f\n", i, fabs(str_1_current - str_1_end), fabs(str_2_current - str_2_end));
#endif

      starting_index = i++;
      // We continue stepping forward in time in increments of 1 while we're still within the bounds of the timespan and either of the two populations
      // are still not within delta of their final occupancy probability. This constitutes the forward step.
    } while (i < num_points - parameters.window_size && (fabs(str_1_current - str_1_end) > delta || fabs(str_2_current - str_2_end) > delta));

    // Since the forward step looks for a position within delta of the *last* user-requested timepoint (to ensure we escape the "warming up" phase), it
    // is posible that the system has already been in equilibrium for some time. The backward step moves us in reverse to find the first position that
    // is not in equilbrium within the user-requested window size. This guarantees that we find the first position that satisfies equilibrium, and not
    // just a position that satisfies (as could happen before, since being within delta of the final position isn't necessarily the first time that you've
    // been in equilibrium; in other words sometimes the forward step moves too far to the right to ensure that we're outside of any local equilibrium).
    while (i >= 0 && is_index_in_equilibrium_within_window_position(eigensystem, parameters, parameters.end_index, i--)) {};

    if (i >= 0 && starting_index != i) {
      if (parameters.verbose) {
        
      }
      printf("Scanned backwards successfully from %d to %d.\n", starting_index, i);

      starting_index = i;
    }
  } else {
    starting_index = -1;
  }

#ifdef DEBUG
  printf("starting_index:\t%d\n", starting_index);
  printf("num_points:\t%d\n", num_points);
#endif

  return starting_index;
}

int is_index_in_equilibrium_within_window_position(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters, int eigensystem_index, int window_start) {
  int i;
  double current_proportion, future_proportion, epsilon = parameters.equilibrium;

  current_proportion = probability_at_time(
                         eigensystem,
                         parameters,
                         pow(10, parameters.start_time + parameters.step_size * window_start),
                         parameters.start_index,
                         eigensystem_index
                       );

  for (i = 1; i < parameters.window_size; ++i) {
    future_proportion = probability_at_time(
                          eigensystem,
                          parameters,
                          pow(10, parameters.start_time + parameters.step_size * (window_start + i)),
                          parameters.start_index,
                          eigensystem_index
                        );

#ifdef INSANE_DEBUG
    printf("%d\t%d\t%d\t%f\t%f\t%e\t%d\n", window_start, window_start + i, eigensystem_index, current_proportion, future_proportion, fabs(current_proportion - future_proportion), fabs(current_proportion - future_proportion) < epsilon);
#endif

    if (fabs(current_proportion - future_proportion) > epsilon) {
#ifdef INSANE_DEBUG
      printf("\n");
#endif

      return 0;
    }
  }

  return 1;
}

void print_equilibrium(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  long double equilibrium_time;

  equilibrium_time = estimate_equilibrium(eigensystem, parameters);

  switch ((int)round(equilibrium_time)) {
    case NEG_INF:
      printf("-Infinity\n");
      break;
    case POS_INF:
      printf("Infinity\n");
      break;
    default:
      printf("%f\n", log(equilibrium_time) / log(10.));
  }
}

void print_population_proportion(const EIGENSYSTEM eigensystem, const POPULATION_PARAMS parameters) {
  double step_counter;

  for (step_counter = parameters.start_time; step_counter <= parameters.end_time; step_counter += parameters.step_size) {
    printf(
      "%+f\t%+.8f\t%+.8f\n",
      step_counter,
      probability_at_time(eigensystem, parameters, pow(10, step_counter), parameters.start_index, parameters.end_index),
      probability_at_time(eigensystem, parameters, pow(10, step_counter), parameters.start_index, parameters.start_index)
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
