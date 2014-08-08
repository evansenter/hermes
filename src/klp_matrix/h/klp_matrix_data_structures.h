#ifndef KLP_MATRIX_DATA_STRUCTURES_H
#define KLP_MATRIX_DATA_STRUCTURES_H

typedef struct {
  int* k;
  int* l;
  double* p;
  int length;
} KLP_MATRIX;

typedef struct {
  double* matrix;
  int row_length;
  char type;
} TRANSITION_MATRIX;

typedef struct {
  int start_state;
  int end_state;
  int max_dist;
  int bp_dist;
  char run_type;
  double epsilon;
  short energy_based;
  short hastings;
  short rate_matrix;
} KLP_PARAMS;

#define TRANSITION_INPUT_FLAG 'T'
#define DIAG_MOVES_ONLY_FLAG 'X'
#define FULLY_CONNECTED_FLAG 'F'
#define RUN_TYPE(run_type, flag) (run_type == flag)
#define MATRIX_TYPE(rate_matrix) (rate_matrix ? 'R' : 'P')

#define T_ROW_ORDER(m, i, j) ((m.matrix)[((i) * (m.row_length)) + (j)])
#define T_COL_ORDER(m, i, j) ((m.matrix)[((j) * (m.row_length)) + (i)])

#endif
