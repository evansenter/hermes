// gcc -lgsl -lgslcblas -M test.c

#include <gsl/gsl_eigen.h>

int main(int argc, char* argv) {
  gsl_eigen_nonsymmv_workspace* workspace = gsl_eigen_nonsymmv_alloc(100);

  return 0;
}
