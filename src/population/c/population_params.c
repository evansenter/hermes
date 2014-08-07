#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "population_constants.h"
#include "population_params.h"

POPULATION_PARAMS init_population_params() {
  POPULATION_PARAMS parameters = {
    .verbose           = 0,
    .input_file        = NULL,
    .sequence          = NULL,
    .start_structure   = NULL,
    .end_structure     = NULL,
    .filename          = NULL,
    .start_index       = -1,
    .end_index         = -1,
    .serialize         = 0,
    .temperature       = 37.,
    .start_time        = -10,
    .end_time          = 10,
    .step_size         = 1e-3,
    .equilibrium       = 0,
    .soft_bounds       = 1,
    .epsilon           = 1e-4,
    .delta             = 1e-3,
    .window_size       = 5,
    .all_subpop_for_eq = 0,
    .lonely_bp         = 0,
    .energy_cap        = 0,
    .eigen_only        = 0,
    .input             = 1,
    .benchmark         = 0
  };
  return parameters;
}

void parse_population_args(POPULATION_PARAMS* parameters, int argc, char** argv, void (*usage)()) {
  int c;
  
  while ((c = getopt(argc, argv, "ogbhnvqa:z:c:s:k:l:i:j:p:e::d:w:t:m:r:f:")) != -1) {
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
        if (!sscanf(optarg, "%d", &parameters->start_index)) {
          (*usage)();
        }
        
        break;
        
      case 'z':
        if (!sscanf(optarg, "%d", &parameters->end_index)) {
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
        if (!sscanf(optarg, "%d", &parameters->serialize)) {
          (*usage)();
        } else if ((int)abs(parameters->serialize) != 1) {
          (*usage)();
        }
        
        break;
        
      case 'f':
        parameters->filename = strdup(optarg);
        break;
        
      case '?':
        switch (optopt) {
          case 'a':
          case 'z':
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
            break;
            
          default:
            if (isprint(optopt)) {
              fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            } else {
              fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            }
        }
        
        (*usage)();
        
      default:
        (*usage)();
    }
  }
  
  if (parameters->verbose) {
    debug_population_parameters(*parameters);
  }
  
  if (population_error_handling(*parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

int population_error_handling(const POPULATION_PARAMS parameters) {
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
  printf("RNAeq parameters\n");
  printf("(c) input_file\t\t\t%s\n",            parameters.input_file);
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
  printf("(m) energy_cap\t\t\t%.1f kcal/mol\n", parameters.energy_cap ? parameters.energy_cap : 10000);
  printf("(q) equilibrium\t\t\t%s\n",           parameters.equilibrium ? "Yes" : "No");
  printf("(e) epsilon\t\t\t%.2e\n",             parameters.epsilon);
  printf("(d) delta\t\t\t%.2e\n",               parameters.delta);
  printf("(n) soft bounds\t\t\t%s\n",           parameters.soft_bounds ? "Yes" : "No");
  printf("(h) all subpop.\t\t\t%s\n",           parameters.all_subpop_for_eq ? "Yes" : "No");
  printf("(w) window size\t\t\t%d\n",           parameters.window_size);
  printf("(g) eigen_only\t\t\t%s\n",            parameters.eigen_only ? "Yes" : "No");
  printf("(b) benchmark\t\t\t%s\n",             parameters.benchmark ? "Yes" : "No");
  
  printf("\n");
}

void population_usage() {
  fprintf(stderr, "RNAeq [options] -c input_csv\n\n");
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");
  fprintf(stderr, "RNAeq [options] -s sequence\n\n");
  fprintf(stderr, "Options include the following:\n");
  population_flags();
  abort();
}

void population_flags() {
  fprintf(stderr, "-a\tstart state,           default is inferred. If provided, should indicate the 0-indexed position in the transition matrix corresponding to the starting structure (see options for -k).\n");
  fprintf(stderr, "-b\t(b)enchmarking,        default is disabled. When enabled, RNAeq will print benchmarking times for internal function calls.\n");
  fprintf(stderr, "-c\t(C)SV input file,      this option is made available to abstain from providing the input CSV as the last command line argument.\n");
  fprintf(stderr, "-d\t(d)elta range,         default is disabled. When provided, this value specifies the delta-size required for the population to be approaching equilibrium. This position is used as a starting point for a more fine-grained scan using the -e and -w values.\n");
  fprintf(stderr, "-e\t(e)psilon,             default is disabled. When provided, we will output the time at which the equilibrium stable within -e for a window of size -w.\n");
  fprintf(stderr, "-f\tbin (f)ilename,        Provided in conjunction with the -r flag to specify the read / write file for serializing the eigensystem.\n");
  fprintf(stderr, "-g\tei(g)envalues only,    default is disabled. When enabled, RNAeq will only generate the eigenvalues for the transition rate matrix.\n");
  fprintf(stderr, "-h\tall subpop. in eq.,    default is disabled. When provided alongside the -q flag, the equilibrium time is estimated for when *all* subpopulations approach equilibrium, rather than when just the start-to-end state curve (-a to -z) doesn't deviate more than epsilon (-e) within the window size (-w).\n");
  fprintf(stderr, "-i\tstart time,            natural log of the starting time for computing population proportion.\n");
  fprintf(stderr, "-j\tend time,              natural log of the ending time for computing population proportion.\n");
  fprintf(stderr, "-k\tstarting structure,    structure for which the probability at time 0 is equal to 1. If not provided, the empty structure is used / inferred.\n");
  fprintf(stderr, "-l\tending structure,      structure of interest for detailed population proportion values. If not provided, the MFE structure is used / inferred.\n");
  fprintf(stderr, "-m\tenergy cap,            default is disabled. When provided, RNAsubopt will only sample structures within energy_cap kcal/mol of the MFE structure.\n");
  fprintf(stderr, "-n\t(n)o soft bounds,      default is disabled. When enabled, the population proportion / equilibrium time will not use soft bounds computed with the -d delta value, and instead will compute population proportion / equilibrium over the entire timespan (timespan can be adjusted with -i and -j).\n");
  fprintf(stderr, "-o\tl(o)nely basepairs,    default is disabled. When enabled, RNAsubopt will sample structures containing lonely base pairs.\n");
  fprintf(stderr, "-p\tste(p) size,           natrual log of the step size for computing population proportion (start_time < start_time + step_size <= end_time).\n");
  fprintf(stderr, "-q\te(q)uilibrium,         Compute the equilibrium time for the target state. Uses an approximation technique that looks at occupancy deviation within a sliding window, parameterized by -e and -w.\n");
  fprintf(stderr, "-r\tse(r)ialize direction, default is disabled. If passed with the value '1', then the eigensystem and its inversion will be serialized to the file specified by the -f flag. If passed with the value '-1', the file specified by the -f flag will be deserialized and used for computing the data of interest.\n");
  fprintf(stderr, "-s\t(s)equence,            sequence of interest for computing population proportions.\n");
  fprintf(stderr, "-t\t(t)emperature,         temperature at which suboptimal structures are generated. This value is passed to (and only used by) ViennaRNA's RNAsubopt.\n");
  fprintf(stderr, "-v\t(v)erbose,             default is disabled, presents some debug information at runtime.\n\n");
  fprintf(stderr, "-w\t(w)indow size,         default is 5. Specifies the window size (exclusive) for predicting equilibrium. Equilibrium is considered as having been achieved when all indices (i + 1)..(i + window_size - 1) are within epsilon (-e) of the population proportion at time i.\n");
  fprintf(stderr, "-z\tend state,             default is inferred. If provided, should indicate the 0-indexed position in the transition matrix corresponding to the ending structure (see options for -l).\n\n");
}
