#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "population_constants.h"
#include "population_params.h"
#include "shared/libklp_matrix_header.h"
#include "shared/constants.h"

POPULATION_PARAMS init_population_params() {
  POPULATION_PARAMS parameters = {
    .input_file         = NULL,
    .sequence           = NULL,
    .start_structure    = NULL,
    .end_structure      = NULL,
    .filename           = NULL,
    .temperature        = 37.,
    .start_time         = -10,
    .end_time           = 10,
    .step_size          = 1e-3,
    .serialize          = 0,
    .equilibrium        = 0,
    .soft_bounds        = 1,
    .epsilon            = 1e-4,
    .delta              = 1e-3,
    .window_size        = 5,
    .num_subpop_to_show = -1,
    .all_subpop_for_eq  = 0,
    .lonely_bp          = 0,
    .energy_cap         = 0,
    .eigen_only         = 0,
    .input              = 1,
    .benchmark          = 0,
    .verbose            = 0
  };
  return parameters;
}

void parse_population_args(KLP_PARAMS* klp_params, POPULATION_PARAMS* parameters, int argc, char** argv, void (*usage)()) {
  int c;
  opterr = 0;
  
  while (optind < argc) {
    if ((c = getopt(argc, argv, "+ogbhnvqa:c:s:k:l:i:j:p:e:d:w:t:m:r:f:")) != -1) {
#ifdef INPUT_DEBUG
      printf("parse_mfpt_args: %c\n", c);
#endif
      
      switch (c) {
        case 'o':
          parameters->lonely_bp = 1;
          break;
          
        case 'g':
          parameters->eigen_only = 1;
          break;
          
        case 'b':
          parameters->benchmark = 1;
          break;
          
        case 'h':
          parameters->all_subpop_for_eq = 1;
          break;
          
        case 'n':
          parameters->soft_bounds = 0;
          break;
          
        case 'v':
          parameters->verbose = 1;
          break;
          
        case 'q':
          parameters->equilibrium = 1;
          break;
          
        case 'c':
          parameters->input_file = strdup(optarg);
          break;
          
        case 's':
          parameters->sequence = strdup(optarg);
          break;
          
        case 'k':
          parameters->start_structure = strdup(optarg);
          break;
          
        case 'l':
          parameters->end_structure = strdup(optarg);
          break;
          
        case 'a':
          if (!sscanf(optarg, "%d", &parameters->num_subpop_to_show)) {
            (*usage)();
          }
          
          break;
          
        case 'i':
          if (!sscanf(optarg, "%lf", &parameters->start_time)) {
            (*usage)();
          }
          
          break;
          
        case 'j':
          if (!sscanf(optarg, "%lf", &parameters->end_time)) {
            (*usage)();
          }
          
          break;
          
        case 'p':
          if (!sscanf(optarg, "%lf", &parameters->step_size)) {
            (*usage)();
          }
          
          break;
          
        case 'e':
          if (!sscanf(optarg, "%lf", &parameters->epsilon)) {
            (*usage)();
          }
          
          parameters->equilibrium = 1;
          break;
          
        case 'd':
          if (!sscanf(optarg, "%lf", &parameters->delta)) {
            (*usage)();
          }
          
          break;
          
        case 'w':
          if (!sscanf(optarg, "%d", &parameters->window_size)) {
            (*usage)();
          }
          
          break;
          
        case 't':
          if (!sscanf(optarg, "%lf", &parameters->temperature)) {
            (*usage)();
          }
          
          break;
          
        case 'm':
          if (!sscanf(optarg, "%lf", &parameters->energy_cap)) {
            (*usage)();
          }
          
          break;
          
        case 'r':
          if (!sscanf(optarg, "%hd", &parameters->serialize)) {
            (*usage)();
          } else if ((int)abs(parameters->serialize) != 1) {
            (*usage)();
          }
          
          break;
          
        case 'f':
          parameters->filename = strdup(optarg);
          break;
          
        case '?':
#ifdef INPUT_DEBUG
          printf("\tcase '?' with %c\n", optopt);
#endif
          
          switch (optopt) {
            case 'a':
            case 's':
            case 'k':
            case 'l':
            case 'e':
            case 'd':
            case 'w':
            case 'i':
            case 'j':
            case 'p':
            case 't':
            case 'c':
            case 'r':
            case 'f':
              fprintf(stderr, "Option -%c requires an argument.\n", optopt);
              (*usage)();
          }
          
          break;
          
        default:
          (*usage)();
      }
    } else {
      optind++;
    }
  }
  
#ifdef INPUT_DEBUG
  printf("Done parsing.\n\n");
#endif
  
  if (parameters->verbose) {
    debug_klp_matrix_parameters(*klp_params);
    debug_population_parameters(*parameters);
  }
  
  if (population_error_handling(*klp_params, *parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

int population_error_handling(const KLP_PARAMS klp_params, const POPULATION_PARAMS parameters) {
  int error = 0;
  
  if (parameters.input) {
    if (((parameters.input_file == NULL) ^ (parameters.sequence == NULL)) != 1) {
      fprintf(stderr, "Error: Either an input file must be specified with -c or an input sequence with -s.\n");
      error++;
    }
    
    if (parameters.input_file != NULL && access(parameters.input_file, R_OK) == -1) {
      fprintf(stderr, "Error: can't read from file %s.\n", parameters.input_file);
      error++;
    }
    
    if (parameters.sequence != NULL) {
      if (parameters.start_structure != NULL && strlen(parameters.sequence) != strlen(parameters.start_structure)) {
        fprintf(stderr, "Error: the starting structure is not the same length as the provided sequence.\n");
        error++;
      }
      
      if (parameters.end_structure != NULL && strlen(parameters.sequence) != strlen(parameters.end_structure)) {
        fprintf(stderr, "Error: the ending structure is not the same length as the provided sequence.\n");
        error++;
      }
    }
  }
  
  if (parameters.num_subpop_to_show < -1) {
    fprintf(stderr, "Error: the number of subpopulations to show must be > 0 (for the top-k subpopulations) or 0 (for all subpopulations).\n");
    error++;
  }
  
  if (parameters.equilibrium && parameters.num_subpop_to_show >= 0 && parameters.all_subpop_for_eq) {
    fprintf(stderr, "Error: requesting to show equilibrium time for multiple subpopulations while forcing all subpopulations to simultaneously be in equilibrium doesn't make sense.\n");
    error++;
  }
  
  if (parameters.energy_cap < 0) {
    fprintf(stderr, "Error: the energy_cap must be a positive number (in kcal/mol) for the energy range above the MFE to sample structures from.\n");
    error++;
  }
  
  if (parameters.epsilon < 0 || parameters.epsilon > 1) {
    fprintf(stderr, "Error: the equilibrium epsilon (provided as the value for the -e flag) must be the epsilon value for considering a population curve to be stable.\n");
    error++;
  }
  
  if (parameters.delta < 0 || parameters.delta > 1) {
    fprintf(stderr, "Error: the equilibrium delta (provided as the value for the -d flag) must be the delta value for considering a population curve to be approaching stability.\n");
    error++;
  }
  
  if (parameters.step_size < 0 || parameters.step_size > 1) {
    fprintf(stderr, "Error: the step size must be between 0 and 1.\n");
    error++;
  }
  
  if (parameters.serialize == -1 && (klp_params.start_state < 0 || klp_params.end_state < 0)) {
    fprintf(stderr, "Error: if you're deserializing, the start and end indices (-A, -Z) in the serialized transition matrix must be provided.\n");
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
  
  if (parameters.window_size < 2) {
    fprintf(stderr, "Error: the window size is exclusive, and must be larger than 2.\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_population_parameters(const POPULATION_PARAMS parameters) {
  printf("RNAeq parameters:\n");
  printf("(c) input_file\t\t\t%s\n", parameters.input_file != NULL ? parameters.input_file : "N/A");
  
  if (parameters.sequence != NULL) {
    printf("(s) sequence\t\t\t%s\n",      parameters.sequence);
    printf("(k) start_structure\t\t%s\n", parameters.sequence != NULL && parameters.start_structure == NULL ? "empty" : parameters.start_structure);
    printf("(l) end_structure\t\t%s\n",   parameters.sequence != NULL && parameters.end_structure == NULL ? "mfe" : parameters.end_structure);
    printf("(t) temperature\t\t\t%.1f\n", parameters.temperature);
    printf("(o) lonely_bp\t\t\t%s\n",     parameters.lonely_bp ? "Yes" : "No");
    
    if (parameters.energy_cap > 0) {
      printf("(m) energy_cap\t\t\t%.1f kcal/mol\n", parameters.energy_cap ? parameters.energy_cap : 10000);
    } else {
      printf("(m) energy_cap\t\t\tN/A\n");
    }
    
  } else {
    printf("(s) sequence\t\t\t%s\n",      "N/A");
    printf("(k) start_structure\t\t%s\n", "N/A");
    printf("(l) end_structure\t\t%s\n",   "N/A");
    printf("(t) temperature\t\t\t%s\n",   "N/A");
    printf("(o) lonely_bp\t\t\t%s\n",     "N/A");
    printf("(m) energy_cap\t\t\t%s\n",    "N/A");
  }
  
  if (parameters.filename != NULL) {
    printf("(f) filename\t\t\t%s\n",  parameters.filename);
    printf("(r) serialize\t\t\t%s\n", (parameters.serialize == -1 ? "writing" : (parameters.serialize == 1 ? "reading" : "unset")));
  } else {
    printf("(f) filename\t\t\t%s\n",  "N/A");
    printf("(r) serialize\t\t\t%s\n", "N/A");
  }
  
  printf("(i) start_time\t\t\t%.2e\n",      parameters.start_time);
  printf("(j) end_time\t\t\t%.2e\n",        parameters.end_time);
  printf("(p) step_size\t\t\t%.2e\n",       parameters.step_size);
  printf("(q) equilibrium\t\t\t%s\n",       parameters.equilibrium ? "Yes" : "No");
  printf("(e) epsilon\t\t\t%.2e\n",         parameters.epsilon);
  printf("(d) delta\t\t\t%.2e\n",           parameters.delta);
  printf("(n) soft bounds\t\t\t%s\n",       parameters.soft_bounds ? "Yes" : "No");
  
  if (parameters.num_subpop_to_show != -1) {
    if (parameters.num_subpop_to_show) {
      printf("(a) all subpop. output\t\ttop %d\n", parameters.num_subpop_to_show);
    } else {
      printf("(a) all subpop. output\t\t%s\n", "all");
    }
  } else {
    printf("(a) all subpop. output\t\t%s\n", "start / end state only");
  }
  
  printf("(h) all subpop. for eq.\t\t%s\n", parameters.all_subpop_for_eq ? "Yes" : "No");
  printf("(w) window size\t\t\t%d\n",       parameters.window_size);
  printf("(g) eigen_only\t\t\t%s\n",        parameters.eigen_only ? "Yes" : "No");
  printf("(b) benchmark\t\t\t%s\n",         parameters.benchmark ? "Yes" : "No");
  
  printf("\n");
}

void population_usage() {
  fprintf(stderr, "RNAeq [options] -s sequence\n");
  fprintf(stderr, "RNAeq [options] -c input_csv\n\n");
  fprintf(stderr, "...where input_csv is a CSV file (with *no* header) of the format:\n\n");
  fprintf(stderr, "k_0,l_0,p(k_0,l_0)\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p(k_n,l_n)\n\n");
  fprintf(stderr, "...with k representing row indices, l representing column indices and p representing the corresponding value at position (k, l) in a row-ordered matrix.\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "------------------------------\n\n");
  klp_matrix_flags();
  population_flags();
  abort();
}

void population_flags() {
  fprintf(stderr, "\t-a\t(a)ll subpop. output,       default is disabled. When enabled, the subsequent option argument indicates number of sub-populations to include in the output, sorted by decreasing popluation occupancy. If the option '0' is provided, all subpopulations are output.\n");
  fprintf(stderr, "\t-b\t(b)enchmarking,             default is disabled. When enabled, will print benchmarking times for internal function calls.\n");
  fprintf(stderr, "\t-c\t(c)SV input file,           this option is made available to abstain from providing the input CSV as the last command line argument.\n");
  fprintf(stderr, "\t-d\t(d)elta range,              default is disabled. When provided, this value specifies the delta-size required for the population to be approaching equilibrium. This position is used as a starting point for a more fine-grained scan using the -e and -w values.\n");
  fprintf(stderr, "\t-e\t(e)psilon,                  default is disabled. When provided, we will output the time at which the equilibrium stable within -e for a window of size -w.\n");
  fprintf(stderr, "\t-f\tbin (f)ilename,             Provided in conjunction with the -r flag to specify the read / write file for serializing the eigensystem.\n");
  fprintf(stderr, "\t-g\tei(g)envalues only,         default is disabled. When enabled, RNAeq will only generate the eigenvalues for the transition rate matrix.\n");
  fprintf(stderr, "\t-h\tall metastates in eq.,      default is disabled. When provided alongside the -q flag, the equilibrium time is estimated for when *all* metastates approach equilibrium, rather than when just the start-to-end state curve (-a to -z) doesn't deviate more than epsilon (-e) within the window size (-w).\n");
  fprintf(stderr, "\t-i\tstart time,                 natural log of the starting time for computing population proportion.\n");
  fprintf(stderr, "\t-j\tend time,                   natural log of the ending time for computing population proportion.\n");
  fprintf(stderr, "\t-k\tstarting structure,         structure for which the probability at time 0 is equal to 1. If not provided, the empty structure is used / inferred.\n");
  fprintf(stderr, "\t-l\tending structure,           structure of interest for detailed population proportion values. If not provided, the MFE structure is used / inferred.\n");
  fprintf(stderr, "\t-m\tenergy cap,                 default is disabled. When provided, RNAsubopt will only sample structures within energy_cap kcal/mol of the MFE structure.\n");
  fprintf(stderr, "\t-n\t(n)o soft bounds,           default is disabled. When enabled, the population proportion / equilibrium time will not use soft bounds computed with the -d delta value, and instead will compute population proportion / equilibrium over the entire timespan (timespan can be adjusted with -i and -j).\n");
  fprintf(stderr, "\t-o\tl(o)nely basepairs,         default is disabled. When enabled, RNAsubopt will sample structures containing lonely base pairs.\n");
  fprintf(stderr, "\t-p\tste(p) size,                natrual log of the step size for computing population proportion (start_time < start_time + step_size <= end_time).\n");
  fprintf(stderr, "\t-q\te(q)uilibrium,              Compute the equilibrium time for the target state. Uses an approximation technique that looks at occupancy deviation within a sliding window, parameterized by -e and -w.\n");
  fprintf(stderr, "\t-r\tse(r)ialize direction,      default is disabled. If passed with the value '1', then the eigensystem and its inversion will be serialized to the file specified by the -f flag. If passed with the value '-1', the file specified by the -f flag will be deserialized and used for computing the data of interest.\n");
  fprintf(stderr, "\t-s\t(s)equence,                 sequence of interest for computing population proportions.\n");
  fprintf(stderr, "\t-t\t(t)emperature,              temperature at which suboptimal structures are generated. This value is passed to (and only used by) ViennaRNA's RNAsubopt.\n");
  fprintf(stderr, "\t-v\t(v)erbose,                  default is disabled, presents some debug information at runtime.\n\n");
  fprintf(stderr, "\t-w\t(w)indow size,              default is 5. Specifies the window size (exclusive) for predicting equilibrium. Equilibrium is considered as having been achieved when all indices (i + 1)..(i + window_size - 1) are within epsilon (-e) of the population proportion at time i.\n\n");
}
