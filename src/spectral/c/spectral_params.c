#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "params.h"

extern double temperature;

SPECTRAL_PARAMS init_spectral_params() {
  SPECTRAL_PARAMS parameters = {
    .verbose          = 0,
    .sequence         = NULL,
    .start_structure  = NULL,
    .end_structure    = NULL,
    .start_index      = -1,
    .end_index        = -1,
    .temperature      = 37.,
    .start_time       = -11,
    .end_time         = 0,
    .step_size        = 1e-1,
    .lonely_bp        = 0,
    .energy_cap       = 0,
    .eigen_only       = 0,
    .benchmark        = 0
  };
  return parameters;
}

SPECTRAL_PARAMS parse_spectral_args(int argc, char** argv) {
  int c;
  SPECTRAL_PARAMS parameters;
  parameters = init_spectral_params();

  while ((c = getopt(argc, argv, "OoGgBbVvS:s:K:k:L:l:I:i:J:j:P:p:T:t:C:c:")) != -1) {
    switch (c) {
      case 'O':
      case 'o':
        parameters.lonely_bp = 1;
        break;

      case 'G':
      case 'g':
        parameters.eigen_only = 1;
        break;

      case 'B':
      case 'b':
        parameters.benchmark = 1;
        break;

      case 'V':
      case 'v':
        parameters.verbose = 1;
        break;

      case 'S':
      case 's':
        parameters.sequence = strdup(optarg);
        break;

      case 'K':
      case 'k':
        parameters.start_structure = strdup(optarg);
        break;

      case 'L':
      case 'l':
        parameters.end_structure = strdup(optarg);
        break;

      case 'I':
      case 'i':
        if (!sscanf(optarg, "%lf", &parameters.start_time)) {
          spectral_usage();
        }

        break;

      case 'J':
      case 'j':
        if (!sscanf(optarg, "%lf", &parameters.end_time)) {
          spectral_usage();
        }

        break;

      case 'P':
      case 'p':
        if (!sscanf(optarg, "%lf", &parameters.step_size)) {
          spectral_usage();
        }

        break;

      case 'T':
      case 't':
        if (!sscanf(optarg, "%lf", &parameters.temperature)) {
          spectral_usage();
        }

        temperature = parameters.temperature;
        break;

      case 'C':
      case 'c':
        if (!sscanf(optarg, "%lf", &parameters.energy_cap)) {
          spectral_usage();
        }

        break;

      case '?':
        switch (optopt) {
          case 'S':
          case 's':
          case 'K':
          case 'k':
          case 'L':
          case 'l':
          case 'E':
          case 'e':
          case 'I':
          case 'i':
          case 'J':
          case 'j':
          case 'P':
          case 'p':
          case 'T':
          case 't':
          case 'C':
          case 'c':
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            break;

          default:
            if (isprint(optopt)) {
              fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            } else {
              fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            }
        }

        spectral_usage();

      default:
        spectral_usage();
    }
  }

  if (parameters.verbose) {
    debug_spectral_parameters(parameters);
  }

  if (spectral_error_handling(parameters)) {
    spectral_usage();
  }

  return parameters;
}

int spectral_error_handling(const SPECTRAL_PARAMS parameters) {
  int error = 0;

  if (parameters.sequence != NULL && parameters.start_structure != NULL && strlen(parameters.sequence) != strlen(parameters.start_structure)) {
    fprintf(stderr, "Error: the starting structure is not the same length as the provided sequence.\n");
    error++;
  }

  if (parameters.sequence != NULL && parameters.end_structure != NULL && strlen(parameters.sequence) != strlen(parameters.end_structure)) {
    fprintf(stderr, "Error: the ending structure is not the same length as the provided sequence.\n");
    error++;
  }

  if (parameters.energy_cap < 0) {
    fprintf(stderr, "Error: the energy_cap must be a positive number (in kcal/mol) for the energy range above the MFE to sample structures from.\n");
    error++;
  }

  if (error) {
    fprintf(stderr, "\n");
  }

  return error;
}

void debug_spectral_parameters(const SPECTRAL_PARAMS parameters) {
  printf("(s) sequence\t\t\t%s\n",              parameters.sequence);
  printf("(k) start_structure\t\t%s\n",         parameters.start_structure == NULL ? "empty" : parameters.start_structure);
  printf("(l) end_structure\t\t%s\n",           parameters.end_structure == NULL ? "mfe" : parameters.end_structure);
  printf("    start_index\t\t\t%d\n",           parameters.start_index);
  printf("    end_index\t\t\t%d\n",             parameters.end_index);
  printf("(t) temperature\t\t\t%.1f\n",         parameters.temperature);
  printf("(i) start_time\t\t\t%.2e\n",          parameters.start_time);
  printf("(j) end_time\t\t\t%.2e\n",            parameters.end_time);
  printf("(p) step_size\t\t\t%.2e\n",           parameters.step_size);
  printf("(o) lonely_bp\t\t\t%s\n",             parameters.lonely_bp ? "Yes" : "No");
  printf("(c) energy_cap\t\t\t%.1f kcal/mol\n", parameters.energy_cap ? parameters.energy_cap : 10000);
  printf("(g) eigen_only\t\t\t%s\n",            parameters.eigen_only ? "Yes" : "No");
  printf("(b) benchmark\t\t\t%s\n",             parameters.benchmark ? "Yes" : "No");
}

void spectral_usage() {
  fprintf(stderr, "RNAspectral [options] -s [sequence]\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-S/s\tsequence (required), the sequence of interest for computing population proportions.\n");
  fprintf(stderr, "-K/k\tstarting structure,  the structure for which the probability at time 0 is equal to 1. If not provided, the empty structure is used.\n");
  fprintf(stderr, "-L/l\tending structure,    the structure of interest for . If not provided, the empty structure is used.\n");
  fprintf(stderr, "-T/t\ttemperature,         the temperature at which suboptimal structures are generated. This value is passed to (and only used by) ViennaRNA's RNAsubopt.\n");
  fprintf(stderr, "-I/i\tstart time,          the natrual log of the starting time for computing population proportion.\n");
  fprintf(stderr, "-J/j\tend time,            the natrual log of the ending time for computing population proportion.\n");
  fprintf(stderr, "-P/p\tstep size,           the natrual log of the step size for computing population proportion (start_time < start_time + step_size <= end_time).\n");
  fprintf(stderr, "-O/o\tlonely basepairs,    the default is disabled. When enabled, RNAsubopt will sample structures containing lonely basepairs .\n");
  fprintf(stderr, "-C/c\tenergy cap,          the default is disabled. When provided, RNAsubopt will only sample structures within energy_cap kcal/mol of the MFE structure.\n");
  fprintf(stderr, "-G/g\teigenvalues only,    the default is disabled. When enabled, RNAspectral will only generate the eigenvalues for the transition rate matrix.\n");
  fprintf(stderr, "-B/b\tbenchmarking,        the default is disabled. When enabled, RNAspectral will print benchmarking times for internal function calls.\n");
  abort();
}
