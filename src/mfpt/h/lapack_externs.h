#ifndef LAPACK_EXTERNS_H
#define LAPACK_EXTERNS_H

#ifdef __cplusplus
extern "C" {
#endif

void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
int dgels_(char* t, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info);

#ifdef __cplusplus
}
#endif

#endif
