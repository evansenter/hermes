#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "mfpt_params.h"
#include "shared/constants.h"

MFPT_PARAMS init_mfpt_params() {
  MFPT_PARAMS parameters = {
    .input_file    = NULL,
    .all_mfpt      = 0,
    .input         = 1,
    .verbose       = 0
  };
  return parameters;
}

void parse_mfpt_args(MFPT_PARAMS* parameters, int argc, char** argv, void (*usage)()) {
  int c;
  
  while ((c = getopt(argc, argv, "lvc:")) != -1) {
    #ifdef INPUT_DEBUG
      printf("parse_mfpt_args: %c\n", c);
    #endif
    
    switch (c) {
      case 'l':
        parameters->all_mfpt = 1;
        break;
        
      case 'v':
        parameters->verbose = 1;
        break;
        
      case 'c':
        parameters->input_file = strdup(optarg);
        break;
        
      case '?':
        #ifdef INPUT_DEBUG
          printf("\tcase '?' with %c\n", optopt);
        #endif
      
        switch (optopt) {
          case 'c':
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            (*usage)();
            
          default:
            break;
        }
        
        break;
        
      default:
        (*usage)();
    }
  }
  
  if (parameters->verbose) {
    debug_mfpt_parameters(*parameters);
  }
  
  if (mfpt_error_handling(*parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

int mfpt_error_handling(const MFPT_PARAMS parameters) {
  int error = 0;
  
  if (parameters.input) {
    // Input CSV exists.
    if (parameters.input_file == NULL) {
      fprintf(stderr, "Error: file not provided. Please provide one as the last argument or with the -c flag.\n");
      error++;
    }
    
    // Input CSV is readable check.
    if (access(parameters.input_file, R_OK) == -1) {
      fprintf(stderr, "Error: can't read from file %s.\n", parameters.input_file);
      error++;
    }
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_mfpt_parameters(const MFPT_PARAMS parameters) {
  printf("RNAmfpt parameters\n");
  printf("(c) input_file\t\t\t%s\n", parameters.input_file);
  printf("(l) all_mfpt\t\t\t%s\n",   parameters.all_mfpt ? "Yes" : "No");
  
  printf("\n");
}

void mfpt_usage() {
  fprintf(stderr, "RNAmfpt [options] -c input_csv\n\n");
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");
  fprintf(stderr, "Options include the following:\n");
  klp_matrix_flags();
  mfpt_flags();
  fprintf(stderr, "Program returns -1 (resp. -2) if the start state (resp. end state) probability is 0. -3 is returned if the distance between the two input structures could not be inferred from the input data (usually also means that one of the states has a 0-probability). Otherwise returns the MFPT as predicted by matrix inversion.\n");
  abort();
}

void mfpt_flags() {
  fprintf(stderr, "-c\t(C)SV input file,           this option is made available to abstain from providing the input CSV as the last command line argument.\n");
  fprintf(stderr, "-l\tprint a(l)l MFPT,           if this flag is provided, the program will print the MFPT for every non-end state to hit the end state, and then print the same beginning -> end MFPT that is printed without this flag. Indices for the MFPT correspond to 0-ordered indices in the transition probability matrix.\n");
  fprintf(stderr, "-v\tverbose,                    default is disabled. If this flag is provided, light debug data will be printed. To enable heavy debugging, use the flags in mfpt_constants.h\n");
}
