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
    .energy_grid_file = NULL,
    .start_structure  = NULL,
    .end_structure    = NULL,
    .start_index      = -1,
    .end_index        = -1,
    .temperature      = 37.,
    .start_time       = -11,
    .end_time         = 0,
    .step_size        = 1e-1,
    .lonely_bp        = 0,
    .energy_cap       = 1,
    .eigen_only       = 0,
    .benchmark        = 0
  };
  return parameters;
}

SPECTRAL_PARAMS parse_spectral_args(int argc, char** argv) {
  int c;
  SPECTRAL_PARAMS parameters;
  parameters = init_spectral_params();
  
  while ((c = getopt(argc, argv, "OoCcGgBbVvS:s:K:k:L:l:E:e:I:i:J:j:P:p:T:t:")) != -1) {
    switch (c) {
      case 'O':
      case 'o':
        parameters.lonely_bp = 1;
        break;
        
      case 'C':
      case 'c':
        parameters.energy_cap = 0;
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
        
      case 'E':
      case 'e':
        parameters.energy_grid_file = strdup(optarg);
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
  
  if (parameters.energy_grid_file != NULL) {
    fprintf(stderr, "Error: energy_grid_file not yet implemented.\n");
    error++;
  }
  
  // Also need error handling for the size of the structures, if provided.
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_spectral_parameters(const SPECTRAL_PARAMS parameters) {
  printf("(s) sequence\t\t\t%s\n",           parameters.sequence);
  printf("(k) start_structure\t\t%s\n",      parameters.start_structure == NULL ? "empty" : parameters.start_structure);
  printf("(l) end_structure\t\t%s\n",        parameters.end_structure == NULL ? "mfe" : parameters.end_structure);
  printf("    start_index\t\t\t%d\n",        parameters.start_index);
  printf("    end_index\t\t\t%d\n",          parameters.end_index);
  printf("(e) energy_grid_file (NYI)\t%s\n", parameters.energy_grid_file == NULL ? "none" : parameters.energy_grid_file);
  printf("(t) temperature\t\t\t%.1f\n",      parameters.temperature);
  printf("(i) start_time\t\t\t%.2e\n",       parameters.start_time);
  printf("(j) end_time\t\t\t%.2e\n",         parameters.end_time);
  printf("(p) step_size\t\t\t%.2e\n",        parameters.step_size);
  printf("(l) lonely_bp\t\t\t%s\n",          parameters.lonely_bp ? "No" : "Yes");
  printf("(c) energy_cap\t\t\t%s\n",         parameters.energy_cap ? "Yes" : "No");
  printf("(g) eigen_only\t\t\t%s\n",         parameters.eigen_only ? "Yes" : "No");
  printf("(b) benchmark\t\t\t%s\n",          parameters.benchmark ? "Yes" : "No");
  printf("(t) temperature\t\t\t%.1f\n",      parameters.temperature);
}

void spectral_usage() {
  fprintf(stderr, "RNAspectral [options] -s [sequence]\n");
  abort();
}
