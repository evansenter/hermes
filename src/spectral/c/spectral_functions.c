#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <tpl.h>
#include "shared/libmfpt_header.h"
#include "constants.h"
#include "initializers.h"
#include "functions.h"

void population_proportion_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS parameters, TRANSITION_MATRIX row_transition_matrix) {
  TRANSITION_MATRIX column_transition_matrix;
  EIGENSYSTEM eigensystem;

  column_transition_matrix = transpose_matrix(row_transition_matrix);
  eigensystem              = convert_transition_matrix_to_eigenvectors(column_transition_matrix);
  invert_matrix(&eigensystem);
  print_population_proportion(parameters, eigensystem);
}

TRANSITION_MATRIX convert_structures_to_transition_matrix(const SOLUTION* all_structures, int num_structures) {
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

void print_population_proportion(const SPECTRAL_PARAMS parameters, const EIGENSYSTEM eigensystem) {
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
