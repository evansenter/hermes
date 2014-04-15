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
    .filename         = NULL,
    .start_index      = -1,
    .end_index        = -1,
    .serialize        = 0,
    .temperature      = 37.,
    .start_time       = -11,
    .end_time         = 0,
    .step_size        = 1e-1,
    .equilibrium      = 0,
    .lonely_bp        = 0,
    .energy_cap       = 0,
    .eigen_only       = 0,
    .input            = 1,
    .benchmark        = 0
  };
  return parameters;
}

void parse_spectral_args(SPECTRAL_PARAMS* parameters, int argc, char** argv) {
  int c;
  
  while ((c = getopt(argc, argv, "OoGgBbVvA:a:Z:z:S:s:K:k:L:l:I:i:J:j:P:p:E:e:T:t:C:c:R:r:F:f:")) != -1) {
    switch (c) {
      case 'O':
      case 'o':
        parameters->lonely_bp = 1;
        break;
        
      case 'G':
      case 'g':
        parameters->eigen_only = 1;
        break;
        
      case 'B':
      case 'b':
        parameters->benchmark = 1;
        break;
        
      case 'V':
      case 'v':
        parameters->verbose = 1;
        break;
        
      case 'S':
      case 's':
        parameters->sequence = strdup(optarg);
        break;
        
      case 'K':
      case 'k':
        parameters->start_structure = strdup(optarg);
        break;
        
      case 'L':
      case 'l':
        parameters->end_structure = strdup(optarg);
        break;
        
      case 'A':
      case 'a':
        if (!sscanf(optarg, "%d", &parameters->start_index)) {
          spectral_usage();
        }
        
        break;
        
      case 'Z':
      case 'z':
        if (!sscanf(optarg, "%d", &parameters->end_index)) {
          spectral_usage();
        }
        
        break;
        
      case 'I':
      case 'i':
        if (!sscanf(optarg, "%lf", &parameters->start_time)) {
          spectral_usage();
        }
        
        break;
        
      case 'J':
      case 'j':
        if (!sscanf(optarg, "%lf", &parameters->end_time)) {
          spectral_usage();
        }
        
        break;
        
      case 'P':
      case 'p':
        if (!sscanf(optarg, "%lf", &parameters->step_size)) {
          spectral_usage();
        }
        
        break;
        
      case 'E':
      case 'e':
        if (!sscanf(optarg, "%lf", &parameters->equilibrium)) {
          spectral_usage();
        }
        
        break;
        
      case 'T':
      case 't':
        if (!sscanf(optarg, "%lf", &parameters->temperature)) {
          spectral_usage();
        }
        
        temperature = parameters->temperature;
        break;
        
      case 'C':
      case 'c':
        if (!sscanf(optarg, "%lf", &parameters->energy_cap)) {
          spectral_usage();
        }
        
        break;
        
      case 'R':
      case 'r':
        if (!sscanf(optarg, "%d", &parameters->serialize)) {
          spectral_usage();
        } else if ((int)abs(parameters->serialize) != 1) {
          spectral_usage();
        }
        
        break;
        
      case 'F':
      case 'f':
        parameters->filename = strdup(optarg);
        break;
        
      case '?':
        switch (optopt) {
          case 'A':
          case 'a':
          case 'Z':
          case 'z':
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
          case 'R':
          case 'r':
          case 'F':
          case 'f':
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
  
  if (parameters->verbose) {
    debug_spectral_parameters(*parameters);
  }
  
  if (spectral_error_handling(*parameters)) {
    spectral_usage();
  }
  
  optind = 1;
}

int spectral_error_handling(const SPECTRAL_PARAMS parameters) {
  int error = 0;
  
  if (parameters.input) {
    if (parameters.sequence != NULL && parameters.start_structure != NULL && strlen(parameters.sequence) != strlen(parameters.start_structure)) {
      fprintf(stderr, "Error: the starting structure is not the same length as the provided sequence.\n");
      error++;
    }
    
    if (parameters.sequence != NULL && parameters.end_structure != NULL && strlen(parameters.sequence) != strlen(parameters.end_structure)) {
      fprintf(stderr, "Error: the ending structure is not the same length as the provided sequence.\n");
      error++;
    }
  }
  
  if (parameters.energy_cap < 0) {
    fprintf(stderr, "Error: the energy_cap must be a positive number (in kcal/mol) for the energy range above the MFE to sample structures from.\n");
    error++;
  }
  
  if (parameters.equilibrium < 0 || parameters.equilibrium > 1) {
    fprintf(stderr, "Error: the equilibrium epsilon (provided as the value for the -e flag) must be the epsilon value for considering a population curve to be stable.\n");
    error++;
  }
  
  if (parameters.serialize == -1 && (parameters.start_index < 0 || parameters.end_index < 0)) {
    fprintf(stderr, "Error: if you're deserializing, the start and end indices (-a, -z) in the serialized transition matrix must be provided.\n");
    error++;
  }
  
  if (parameters.serialize == 1 && parameters.eigen_only) {
    fprintf(stderr, "Error: you can't serialize to a file (-r) and compute only the eigenvalues (-g). This is because the serialized data structure has the inverted eigenvector-matrix in it, which is not computed with the -g flag for performance reasons.\n");
    error++;
  }
  
  if ((parameters.serialize && parameters.filename == NULL) || (!parameters.serialize && parameters.filename != NULL)) {
    fprintf(stderr, "Error: can't serialize (-r) without a filename provided (-f).\n");
    error++;
  }
  
  if (parameters.filename && parameters.serialize) {
    if (parameters.serialize == 1 && access(parameters.filename, F_OK) != -1) {
      fprintf(stderr, "Error: file %s exists, cowering out.\n", parameters.filename);
      error++;
    }
    
    if (parameters.serialize == -1 && access(parameters.filename, R_OK) == -1) {
      fprintf(stderr, "Error: can't read from file %s.\n", parameters.filename);
      error++;
    }
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_spectral_parameters(const SPECTRAL_PARAMS parameters) {
  printf("RNAspectral parameters\n");
  printf("(s) sequence\t\t\t%s\n",              parameters.sequence);
  printf("(k) start_structure\t\t%s\n",         parameters.start_structure == NULL && parameters.start_index < 0 ? "empty" : parameters.start_structure);
  printf("(l) end_structure\t\t%s\n",           parameters.end_structure == NULL && parameters.end_index < 0 ? "mfe" : parameters.end_structure);
  printf("(f) filename\t\t\t%s\n",              parameters.filename);
  printf("(a) start_index\t\t\t%d\n",           parameters.start_index);
  printf("(z) end_index\t\t\t%d\n",             parameters.end_index);
  printf("(r) serialize\t\t\t%d\n",             parameters.serialize);
  printf("(t) temperature\t\t\t%.1f\n",         parameters.temperature);
  printf("(i) start_time\t\t\t%.2e\n",          parameters.start_time);
  printf("(j) end_time\t\t\t%.2e\n",            parameters.end_time);
  printf("(p) step_size\t\t\t%.2e\n",           parameters.step_size);
  printf("(o) lonely_bp\t\t\t%s\n",             parameters.lonely_bp ? "Yes" : "No");
  printf("(c) energy_cap\t\t\t%.1f kcal/mol\n", parameters.energy_cap ? parameters.energy_cap : 10000);
  printf("(e) equilibrium\t\t\t%.2e\n",         parameters.equilibrium);
  printf("(g) eigen_only\t\t\t%s\n",            parameters.eigen_only ? "Yes" : "No");
  printf("(b) benchmark\t\t\t%s\n",             parameters.benchmark ? "Yes" : "No");
  
  printf("\n");
}

void spectral_usage() {
  fprintf(stderr, "RNAspectral [options] -s [sequence]\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-A/a\tstart state,           default is inferred. If provided, should indicate the 0-indexed position in the transition matrix corresponding to the starting structure (see options for -k).\n");
  fprintf(stderr, "-B/b\t(b)enchmarking,        default is disabled. When enabled, RNAspectral will print benchmarking times for internal function calls.\n");
  fprintf(stderr, "-C/c\tenergy (c)ap,          default is disabled. When provided, RNAsubopt will only sample structures within energy_cap kcal/mol of the MFE structure.\n");
  fprintf(stderr, "-E/e\t(e)quilibrium,         default is disabled. When provided, we will output the predicted equilibrium time.\n");
  fprintf(stderr, "-G/g\tei(g)envalues only,    default is disabled. When enabled, RNAspectral will only generate the eigenvalues for the transition rate matrix.\n");
  fprintf(stderr, "-I/i\tstart time,            natrual log of the starting time for computing population proportion.\n");
  fprintf(stderr, "-J/j\tend time,              natrual log of the ending time for computing population proportion.\n");
  fprintf(stderr, "-K/k\tstarting structure,    structure for which the probability at time 0 is equal to 1. If not provided, the empty structure is used / inferred.\n");
  fprintf(stderr, "-L/l\tending structure,      structure of interest for detailed population proportion values. If not provided, the MFE structure is used / inferred.\n");
  fprintf(stderr, "-O/o\tl(o)nely basepairs,    default is disabled. When enabled, RNAsubopt will sample structures containing lonely basepairs.\n");
  fprintf(stderr, "-P/p\tste(p) size,           natrual log of the step size for computing population proportion (start_time < start_time + step_size <= end_time).\n");
  fprintf(stderr, "-R/r\tse(r)ialize direction, default is disabled. If passed with the value '1', then the eigensystem and its inversion will be serialized to the file specified by the -f flag. If passed with the value '-1', the file specified by the -f flag will be deserialized and used for computing the data of interest.\n");
  fprintf(stderr, "-S/s\t(s)equence,            sequence of interest for computing population proportions.\n");
  fprintf(stderr, "-T/t\t(t)emperature,         temperature at which suboptimal structures are generated. This value is passed to (and only used by) ViennaRNA's RNAsubopt.\n");
  fprintf(stderr, "-Z/z\tend state,             default is inferred. If provided, should indicate the 0-indexed position in the transition matrix corresponding to the ending structure (see options for -l).\n");
  abort();
}
