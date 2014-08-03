#ifndef KLP_MATRIX_DATA_STRUCTURES_H
#define KLP_MATRIX_DATA_STRUCTURES_H

#define T_ROW_ORDER(m, i, j) ((m.matrix)[((i) * (m.row_length)) + (j)])
#define T_COL_ORDER(m, i, j) ((m.matrix)[((j) * (m.row_length)) + (i)])

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

#endif
