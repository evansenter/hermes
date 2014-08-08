#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <limits>
#include "fftbor2d_params.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_LENGTH 1024
#define TRIEQUALS(x, y, z) ((x == y) && (y == z)) /* Transitivity (BOOM) */

FFTBOR2D_PARAMS init_fftbor2d_params() {
  FFTBOR2D_PARAMS parameters = {
    NULL,                                  // sequence
    NULL,                                  // structure_1
    NULL,                                  // structure_2
    NULL,                                  // energy_file
    0,                                     // seq_length
    (int)ceil(log(pow(10., 8)) / log(2.)), // precision
    0,                                     // max_threads
    BASIC_FLAG,                            // format
    0,                                     // benchmark
    0                                      // verbose
  };
#ifdef _OPENMP
  parameters.max_threads = omp_get_max_threads();
#endif
  return parameters;
}

void free_fftbor2d_params(FFTBOR2D_PARAMS parameters) {
  free(parameters.energy_file);
}

void parse_fftbor2d_args(FFTBOR2D_PARAMS& parameters, int argc, char** argv, void (*usage)()) {
  int c;
  
  while ((c = getopt(argc, argv, "VvBbMmSsCcT:t:E:e:I:i:J:j:K:k:P:p:")) != -1) {
    switch (c) {
      case 'V':
      case 'v':
        parameters.verbose = 1;
        break;
        
      case 'B':
      case 'b':
        parameters.benchmark = 1;
        break;
        
      case 'M':
      case 'm':
        parameters.format = MATRIX_FLAG;
        break;
        
      case 'S':
      case 's':
        parameters.format = SIMPLE_FLAG;
        break;
        
      case 'C':
      case 'c':
        parameters.format = CSV_FLAG;
        break;
        
      case 'T':
      case 't':
        if (!sscanf(optarg, "%lf", &temperature)) {
          (*usage)();
        }
        
        break;
        
      case 'P':
      case 'p':
        if (!sscanf(optarg, "%d", &parameters.precision)) {
          (*usage)();
        } else if (parameters.precision < 0 || parameters.precision > std::numeric_limits<double>::digits) {
          (*usage)();
        }
        
        break;
        
      case 'E':
      case 'e':
        parameters.energy_file = strdup(optarg);
        break;
        
        
      case 'I':
      case 'i':
        parameters.sequence   = strdup(optarg);
        parameters.seq_length = strlen(parameters.sequence);
        break;
        
      case 'J':
      case 'j':
        parameters.structure_1 = strdup(optarg);
        break;
        
      case 'K':
      case 'k':
        parameters.structure_2 = strdup(optarg);
        break;
        
      case '?':
        switch (optopt) {
          case 'T':
          case 't':
          case 'P':
          case 'p':
          case 'E':
          case 'e':
          case 'I':
          case 'i':
          case 'J':
          case 'j':
          case 'K':
          case 'k':
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
  
  if (parameters.sequence == NULL || parameters.structure_1 == NULL || parameters.structure_2 == NULL) {
    if (optind == argc) {
      (*usage)();
    } else {
      parse_fftbor2d_sequence_data(argc, argv, optind, parameters, &fftbor2d_usage);
    }
  }
  
  if (parameters.verbose) {
    debug_fftbor2d_parameters(parameters);
  }
  
  if (fftbor2d_error_handling(parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

void parse_fftbor2d_sequence_data(int argc, char** argv, int argp, FFTBOR2D_PARAMS& parameters, void (*usage)()) {
  int i, char_index = 0;
  FILE* file;
  char line[MAX_LENGTH];
  file = fopen(argv[argp], "r");
  
  if (file == NULL) {
    /* Input is not a file */
    /* argv[argp] should be sequence and argv[argp + 1], argv[argp + 2] should be structures */
    if (argc <= argp + 2) {
      (*usage)();
    }
    
    parameters.seq_length  = strlen(argv[argp]);
    parameters.sequence    = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_1 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_2 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    sscanf(argv[argp++], "%s", parameters.sequence);
    sscanf(argv[argp++], "%s", parameters.structure_1);
    sscanf(argv[argp],   "%s", parameters.structure_2);
  } else {
    /* Input is a file */
    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading file\n");
      exit(0);
    }
    
    while (*line == '*' || *line == '\0' || *line == '>') {
      if (fgets(line, sizeof(line), file) == NULL) {
        break;
      }
    }
    
    if (line == NULL) {
      (*usage)();
    }
    
    // This was a tricky bug to catch, fgets (perhaps obviously) reads a line in, including \n character,
    // which is different from terminating character \0, making data.seq_length off by one.
    while (line[char_index] != '\0' && char_index < MAX_LENGTH) {
      if (line[char_index] == '\n') {
        line[char_index] = '\0';
      }
      
      char_index++;
    }
    
    parameters.seq_length  = strlen(line);
    parameters.sequence    = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_1 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_2 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    sscanf(line, "%s", parameters.sequence);
    
    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading file\n");
      exit(0);
    }
    
    sscanf(line, "%s", parameters.structure_1);
    
    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading file\n");
      exit(0);
    }
    
    sscanf(line, "%s", parameters.structure_2);
    fclose(file);
  }
  
  parameters.sequence[parameters.seq_length]    = '\0';
  parameters.structure_1[parameters.seq_length] = '\0';
  parameters.structure_2[parameters.seq_length] = '\0';
  
  /* Convert RNA sequence to uppercase and make sure there are no Ts in sequence (replace by U). */
  for (i = 0; i < parameters.seq_length; ++i) {
    parameters.sequence[i] = toupper(parameters.sequence[i]);
    
    if (parameters.sequence[i] == 'T') {
      parameters.sequence[i] = 'U';
    }
  }
}

int fftbor2d_error_handling(const FFTBOR2D_PARAMS parameters) {
  int error = 0;
  
  if (!TRIEQUALS(strlen(parameters.sequence), strlen(parameters.structure_1), strlen(parameters.structure_2))) {
    fprintf(
      stderr,
      "(%d) %s\n(%d) %s\n(%d) %s\n",
      (int)strlen(parameters.sequence),
      parameters.sequence,
      (int)strlen(parameters.structure_1),
      parameters.structure_1,
      (int)strlen(parameters.structure_2),
      parameters.structure_2
    );
    fprintf(stderr, "Length of RNA sequence and structures must be equal.\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_fftbor2d_parameters(const FFTBOR2D_PARAMS parameters) {
  printf("FFTbor2D parameters\n");
  printf("(i) sequence\t\t\t%s\n",      parameters.sequence);
  printf("(j) structure_1\t\t\t%s\n",   parameters.structure_1);
  printf("(k) structure_2\t\t\t%s\n",   parameters.structure_2);
  printf("    seq_length\t\t\t%d\n",    parameters.seq_length);
  printf("    max_threads\t\t\t%d\n",   parameters.max_threads);
  printf("    format\t\t\t%c\n",        parameters.format);
  printf("(t) temperature\t\t\t%.1f\n", temperature);
  printf("(p) precision\t\t\t%d\n",     parameters.precision);
  printf("(e) energy_file\t\t\t%s\n",   parameters.energy_file);
  printf("(b) benchmark\t\t\t%s\n",     parameters.benchmark ? "Yes" : "No");
  
  printf("\n");
}

void fftbor2d_usage() {
  fprintf(stderr, "FFTbor2D [options] sequence structure_1 structure_2\n\n");
  fprintf(stderr, "FFTbor2D [options] -i sequence -j structure_1 -k structure_2\n\n");
  fprintf(stderr, "FFTbor2D [options] filename\n");
  fprintf(stderr, "where filename is a file of format:\n");
  fprintf(stderr, "\t>comment (optional line)\n");
  fprintf(stderr, "\tsequence (max length: %d)\n", MAX_LENGTH);
  fprintf(stderr, "\tsecondary structure (1)\n");
  fprintf(stderr, "\tsecondary structure (2)\n\n");
  fprintf(stderr, "Options include the following:\n");
  fftbor2d_flags();
  fprintf(stderr, "Note: output formatting flags (C/c, M/m, S/s) are mutually exclusive. If more than one is provided, *only* last flag will be honored.\n");
  abort();
}

void fftbor2d_flags() {
  fprintf(stderr, "\t-B/b\t(b)enchmark,     default is off. If on, benchmarking data will print alongside normal results.\n");
  fprintf(stderr, "\t-C/c\t(C)SV output,    default is disabled, presents output in CSV format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z).\n");
  fprintf(stderr, "\t-E/e\t(e)nergyfile,    default is rna_turner2004.par in this current directory. Must be name of a file with all energy parameters (in same format as used in Vienna RNA). Energy file lookup first checks current directory, and then iterates through PATH shell variable until a matching file is found. If no file is found, default ViennaRNA parameters are used and a warning is presented to user. If -E switch is explicitly provided, that file is used in lieu of searching for rna_turner2004.par file.\n");
  fprintf(stderr, "\t-I/i\tsequence,        The sequence to be used by FFTbor2D.\n");
  fprintf(stderr, "\t-J/j\tstructure_1,     The first structure to use with FFTbor2D.\n");
  fprintf(stderr, "\t-K/k\tstructure_2,     The second structure to use with FFTbor2D.\n");
  fprintf(stderr, "\t-M/m\t(m)atrix format, default is disabled, presents output in a matrix format instead of a column format.\n");
  fprintf(stderr, "\t-P/p\t(p)recision,     default is %d, indicates precision (base 2) of probabilities Z_k / Z to be returned (0-%d, 0 disables precision handling).\n", (int)ceil(log(pow(10., 8)) / log(2.)), std::numeric_limits<double>::digits);
  fprintf(stderr, "\t-S/s\t(s)imple output, default is disabled, presents output in column format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z), -RTln(Z_{k,l})).\n");
  fprintf(stderr, "\t-T/t\t(t)emperature,   default is 37 degrees Celsius.\n");
  fprintf(stderr, "\t-V/v\t(v)erbose,       default is disabled, presents some debug information at runtime.\n\n");
}
