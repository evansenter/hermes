#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "klp_matrix_params.h"

KLP_PARAMS init_klp_matrix_params() {
  KLP_PARAMS parameters = {
    .start_state   = -1,
    .end_state     = -1,
    .bp_dist       = 0,
    .max_dist      = 0,
    .epsilon       = 0.,
    .run_type      = DIAG_MOVES_ONLY_FLAG,
    .energy_based  = 0,
    .hastings      = 0,
    .rate_matrix   = 0,
    .verbose       = 0
  };
  return parameters;
}

void parse_klp_matrix_args(KLP_PARAMS* parameters, int argc, char** argv, void (*usage)()) {
  int c, option_index = 0;
  
  while (1) {
    static struct option long_options[] = {
      { "start_state",       required_argument, 0, 'A' },
      { "end_state",         required_argument, 0, 'Z' },
      { "bp_dist",           required_argument, 0, 'D' },
      { "max_dist",          required_argument, 0, 'N' },
      { "epsilon",           required_argument, 0, 'O' },
      { "fully_connected",   no_argument,       0, 'F' },
      { "transition_matrix", no_argument,       0, 'T' },
      { "single_bp_only",    no_argument,       0, 'X' },
      { "energy_based",      no_argument,       0, 'E' },
      { "hastings",          no_argument,       0, 'H' },
      { "rate_matrix",       no_argument,       0, 'R' },
      { "verbose",           no_argument,       0, 'V' },
      { 0,                   0,                 0,  0  }
    };
    
    c = getopt_long(argc, argv, "A:Z:D:N:O:FTXEHRV", long_options, &option_index);
    
    if (c == -1) {
      break;
    }
    
    switch (c) {
      case 0:
        break;
        
      case 'F':
        parameters->run_type = FULLY_CONNECTED_FLAG;
        break;
        
      case 'T':
        parameters->run_type = TRANSITION_INPUT_FLAG;
        break;
        
      case 'X':
        parameters->run_type = DIAG_MOVES_ONLY_FLAG;
        break;
        
      case 'E':
        parameters->energy_based = 1;
        break;
        
      case 'H':
        parameters->hastings = 1;
        break;
        
      case 'R':
        parameters->rate_matrix = 1;
        break;
        
      case 'V':
        parameters->verbose = 1;
        break;
        
      case 'A':
        if (!sscanf(optarg, "%d", &parameters->start_state)) {
          (*usage)();
        } else if (parameters->start_state < 0) {
          (*usage)();
        }
        
        break;
        
      case 'Z':
        if (!sscanf(optarg, "%d", &parameters->end_state)) {
          (*usage)();
        } else if (parameters->end_state < 0) {
          (*usage)();
        }
        
        break;
        
      case 'N':
        if (!sscanf(optarg, "%d", &parameters->max_dist)) {
          (*usage)();
        } else if (parameters->max_dist <= 0) {
          (*usage)();
        }
        
        break;
        
      case 'D':
        if (!sscanf(optarg, "%d", &parameters->bp_dist)) {
          (*usage)();
        } else if (parameters->bp_dist <= 0) {
          (*usage)();
        }
        
        break;
        
      case 'O':
        if (!sscanf(optarg, "%lf", &parameters->epsilon)) {
          (*usage)();
        } else if (parameters->epsilon <= 1e-16) {
          (*usage)();
        }
        
        break;
        
      case '?':
        (*usage)();
        
      default:
        (*usage)();
    }
  }
  
  if (parameters->verbose) {
    debug_klp_matrix_parameters(*parameters);
  }
  
  if (klp_matrix_error_handling(*parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

int klp_matrix_error_handling(const KLP_PARAMS parameters) {
  int error = 0;
  
  // Type of run check.
  if (RUN_TYPE(parameters, TRANSITION_INPUT_FLAG) + RUN_TYPE(parameters, DIAG_MOVES_ONLY_FLAG) + RUN_TYPE(parameters, FULLY_CONNECTED_FLAG) != 1) {
    fprintf(stderr, "Error: exactly one of -t, -x or -f must be provided!\n");
    error++;
  }
  
  // Transition matrix input requirements.
  if (RUN_TYPE(parameters, TRANSITION_INPUT_FLAG) && !(parameters.start_state >= 0 && parameters.end_state >= 0)) {
    fprintf(stderr, "Error: if the -t flag is provided, -a and -z must be explicitly set!\n");
    error++;
  }
  
  // Transition matrix input restrictions.
  if (RUN_TYPE(parameters, TRANSITION_INPUT_FLAG) && (parameters.hastings || parameters.energy_based)) {
    fprintf(stderr, "Error: if the -t flag is provided, -h and -e are not permitted!\n");
    error++;
  }
  
  // Fully connected graph restrictions.
  if (RUN_TYPE(parameters, FULLY_CONNECTED_FLAG) && (parameters.hastings || parameters.energy_based)) {
    fprintf(stderr, "Error: if the -f flag is provided, -h and -e are not permitted (or -a and -z without -d)!\n");
    error++;
  }
  
  // Extending k/l/p restrictions.
  if (parameters.max_dist && parameters.energy_based) {
    fprintf(stderr, "Error: if the -n flag is provided, -e is not permitted!\n");
    error++;
  }
  
  // If we're extending k/l/p we need to know how the user wants the zero-positions filled.
  if (parameters.max_dist && !parameters.epsilon) {
    fprintf(stderr, "Error: if using the full grid (bounded by -n), -o needs to be specified for populating 0-probability positions!\n");
    error++;
  }
  
  // If we're willing the zero-positions, we need to know the bounding size.
  if (parameters.epsilon && !parameters.max_dist) {
    fprintf(stderr, "Error: if using -o, the full grid (bounded by -n) needs to be specified for populating 0-probability positions!\n");
    error++;
  }
  
  // MFPT will be 0.
  if (parameters.start_state == parameters.end_state && parameters.start_state >= 0) {
    fprintf(stderr, "Error: if the -a and -z flags are identical the MFPT is 0!\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_klp_matrix_parameters(const KLP_PARAMS parameters) {
  char* buffer = calloc(128, sizeof(char));
  
  printf("CSV parsing parameters\n");
  printf("(F) fully_connected\t\t%s\n",       RUN_TYPE(parameters, FULLY_CONNECTED_FLAG)  ? "Yes" : "No");
  printf("(T) transition_matrix_input\t%s\n", RUN_TYPE(parameters, TRANSITION_INPUT_FLAG) ? "Yes" : "No");
  printf("(X) single_bp_moves_only\t%s\n",    RUN_TYPE(parameters, DIAG_MOVES_ONLY_FLAG)  ? "Yes" : "No");
  printf("(E) energy_based\t\t%s\n",          parameters.energy_based                     ? "Yes" : "No");
  printf("(H) hastings\t\t\t%s\n",            parameters.hastings                         ? "Yes" : "No");
  printf("(R) rate_matrix\t\t\t%s\n",         parameters.rate_matrix                      ? "Yes" : "No");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.start_state);
  printf("(A) start_state\t\t\t%s\n", parameters.start_state >= 0 ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.end_state);
  printf("(Z) end_state\t\t\t%s\n", parameters.end_state >= 0 ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.max_dist);
  printf("(N) max_dist\t\t\t%s\n", parameters.max_dist ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.bp_dist);
  printf("(D) bp_dist\t\t\t%s\n", parameters.bp_dist ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%.2e", parameters.epsilon);
  printf("(O) epsilon\t\t\t%s\n", parameters.epsilon >= 0 ? buffer : "N/A");
  
  printf("\n");
}

void klp_matrix_usage() {
  klp_matrix_flags();
  abort();
}

void klp_matrix_flags() {
  fprintf(stderr, "-A\tstart state,                default is -1 (inferred from input data as the first row in the CSV whose entry in the first column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the start state.\n");
  fprintf(stderr, "-D\tstart/end (d)istance,       default is disabled. When provided, indicates the base pair distance between the starting / ending structures. This flag is used in conjunction with the -n flag, and is needed in cases when the base pair distance between the two structures can't be inferred from the input grid.\n");
  fprintf(stderr, "-E\t(e)nergy-based transitions, default is disabled. If this flag is provided, the transition from state a to b will be calculated as (min(1, exp(-(E_b - E_a) / RT)) / n) rather than (min(1, p_b / p_a) / n).\n");
  fprintf(stderr, "-F\t(f)ully connected,          if the graph is fully connected we permit transitions between arbitrary positions in the 2D grid.\n");
  fprintf(stderr, "-H\t(H)astings adjustment,      default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted (in the all-to-all transition case, N(X) / N(Y) == 1). Calculating N(X) and N(Y) will respect grid boundaries and the triangle equality, and the basepair distance between the two structures for kinetics is inferred from the energy grid.\n");
  fprintf(stderr, "-N\tsequence le(n)gth,          default is disabled. This flag represents the sequence length of the sequence on which kinetics is being performed. It is used to ensure that the graph is fully connected.\n");
  fprintf(stderr, "-O\tepsil(o)n,                  if the graph is going to be populated with all possible moves (via the -n flag), this will inflate all 0-probability positions.\n");
  fprintf(stderr, "-R\t(r)ate matrix,              default is disabled. If this flag is provided, the transition rate matrix is computed rather than the transition probability matrix.\n");
  fprintf(stderr, "-T\t(t)ransition matrix input,  default is disabled. If this flag is provided, the input is expected to be a transition probability matrix, rather than a 2D energy grid. In this case, the first two columns in the CSV file are row-order indices into the transition probability matrix, and the third (final) column is the transition probability of that cell.\n");
  fprintf(stderr, "-V\tverbose,                    default is disabled. If this flag is provided, light debug data will be printed. To enable heavy debugging, use the flags in mfpt_constants.h\n");
  fprintf(stderr, "-X\tsingle basepair moves,      default is enabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted. This option makes the assumption that the input is *not* a transition probability matrix already, and the input energy grid already satisfies the triangle inequality / parity condition.\n");
  fprintf(stderr, "-Z\tend state,                  default is -1 (inferred from input data as the first row in the CSV whose entry in the second column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the end state.\n\n");
}
