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
    NULL,                                  // filename
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
  opterr = 0;
  
  while (optind < argc) {
    if ((c = getopt(argc, argv, "+vbmsct:e:f:i:j:k:p:")) != -1) {
      #ifdef INPUT_DEBUG
        printf("parse_mfpt_args: %c\n", c);
      #endif
      
      switch (c) {
        case 'v':
          parameters.verbose = 1;
          break;
        
        case 'b':
          parameters.benchmark = 1;
          break;
        
        case 'm':
          parameters.format = MATRIX_FLAG;
          break;
        
        case 's':
          parameters.format = SIMPLE_FLAG;
          break;
        
        case 'c':
          parameters.format = CSV_FLAG;
          break;
        
        case 't':
          if (!sscanf(optarg, "%lf", &temperature)) {
            (*usage)();
          }
        
          break;
        
        case 'p':
          if (!sscanf(optarg, "%d", &parameters.precision)) {
            (*usage)();
          } else if (parameters.precision < 0 || parameters.precision > std::numeric_limits<double>::digits) {
            (*usage)();
          }
        
          break;
        
        case 'e':
          parameters.energy_file = strdup(optarg);
          break;
        
        
        case 'f':
          parameters.input_file = strdup(optarg);
          break;
        
        case 'i':
          parameters.sequence   = strdup(optarg);
          parameters.seq_length = strlen(parameters.sequence);
          break;
        
        case 'j':
          parameters.structure_1 = strdup(optarg);
          break;
        
        case 'k':
          parameters.structure_2 = strdup(optarg);
          break;
        
        case '?':
          #ifdef INPUT_DEBUG
            printf("\tcase '?' with %c\n", optopt);
          #endif
  
          switch (optopt) {
            case 't':
            case 'p':
            case 'e':
            case 'i':
            case 'j':
            case 'k':
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
  
  if (parameters.input_file != NULL) {
    parse_fftbor2d_data_from_file(parameters, &fftbor2d_usage);
  }
  
  
  
  if (parameters.verbose) {
    debug_fftbor2d_parameters(parameters);
  }
  
  if (fftbor2d_error_handling(parameters)) {
    (*usage)();
  }
  
  optind = 1;
}

void parse_fftbor2d_data_from_file(FFTBOR2D_PARAMS& parameters, void (*usage)()) {
  int char_index = 0;
  FILE* file;
  char line[MAX_LENGTH];
  
  if (access(parameters.input_file, R_OK) == -1) {
    fprintf(stderr, "Error: can't read from file %s.\n", parameters.input_file);
    (*usage)();
  }
  
  file = fopen(parameters.input_file, "r");
  
  if (fgets(line, sizeof(line), file) == NULL) {
    fprintf(stderr, "Error: There was an error while reading file %s.\n", parameters.input_file);
    (*usage)();
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
    fprintf(stderr, "Error: There was an error while reading file %s.\n", parameters.input_file);
    (*usage)();
  }
  
  sscanf(line, "%s", parameters.structure_1);
  
  if (fgets(line, sizeof(line), file) == NULL) {
    fprintf(stderr, "Error: There was an error while reading file %s.\n", parameters.input_file);
    (*usage)();
  }
  
  sscanf(line, "%s", parameters.structure_2);
  fclose(file);
  
  parameters.sequence[parameters.seq_length]    = '\0';
  parameters.structure_1[parameters.seq_length] = '\0';
  parameters.structure_2[parameters.seq_length] = '\0';
}

void sanitize_sequence_and_structure_parameters(FFTBOR2D_PARAMS& parameters) {
  int i;
  
  for (i = 0; i < parameters.seq_length; ++i) {
    parameters.sequence[i] = toupper(parameters.sequence[i]);
    
    if (parameters.sequence[i] == 'T') {
      parameters.sequence[i] = 'U';
    }
  }
}

int fftbor2d_error_handling(const FFTBOR2D_PARAMS parameters) {
  int error = 0;
  
  if (TRIEQUALS(parameters.input_file, parameters.sequence, NULL)) {
    fprintf(stderr, "Error: either an input file must be specified with -f or the sequence and structures via -i, -j, -k.\n");
    error++;
  } else {
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
      fprintf(stderr, "Error: length of RNA sequence and structures must be equal.\n");
      error++;
    }
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_fftbor2d_parameters(const FFTBOR2D_PARAMS parameters) {
  printf("FFTbor2D parameters:\n");
  printf("(f) filename\t\t\t%s\n",      parameters.input_file != NULL ? parameters.input_file : "N/A");
  printf("(i) sequence\t\t\t%s\n",      parameters.sequence);
  printf("(j) structure_1\t\t\t%s\n",   parameters.structure_1);
  printf("(k) structure_2\t\t\t%s\n",   parameters.structure_2);
  printf("    seq_length\t\t\t%d\n",    parameters.seq_length);
  printf("    max_threads\t\t\t%d\n",   parameters.max_threads);
  printf("    format\t\t\t%s\n",        HUMANIZED_FORMAT);
  printf("(t) temperature\t\t\t%.1f\n", temperature);
  printf("(p) precision\t\t\t%d\n",     parameters.precision);
  printf("(e) energy_file\t\t\t%s\n",   parameters.energy_file != NULL ? parameters.energy_file : "not specified");
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
  fprintf(stderr, "\t-b\t(b)enchmark,     default is off. If on, benchmarking data will print alongside normal results.\n");
  fprintf(stderr, "\t-c\t(C)SV output,    default is disabled, presents output in CSV format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z).\n");
  fprintf(stderr, "\t-e\t(e)nergyfile,    default is rna_turner2004.par in this current directory. Must be name of a file with all energy parameters (in same format as used in Vienna RNA). Energy file lookup first checks current directory, and then iterates through PATH shell variable until a matching file is found. If no file is found, default ViennaRNA parameters are used and a warning is presented to user. If -E switch is explicitly provided, that file is used in lieu of searching for rna_turner2004.par file.\n");
  fprintf(stderr, "\t-i\tsequence,        The sequence to be used by FFTbor2D.\n");
  fprintf(stderr, "\t-j\tstructure_1,     The first structure to use with FFTbor2D.\n");
  fprintf(stderr, "\t-k\tstructure_2,     The second structure to use with FFTbor2D.\n");
  fprintf(stderr, "\t-m\t(m)atrix format, default is disabled, presents output in a matrix format instead of a column format.\n");
  fprintf(stderr, "\t-p\t(p)recision,     default is %d, indicates precision (base 2) of probabilities Z_k / Z to be returned (0-%d, 0 disables precision handling).\n", (int)ceil(log(pow(10., 8)) / log(2.)), std::numeric_limits<double>::digits);
  fprintf(stderr, "\t-s\t(s)imple output, default is disabled, presents output in column format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z), -RTln(Z_{k,l})).\n");
  fprintf(stderr, "\t-t\t(t)emperature,   default is 37 degrees Celsius.\n");
  fprintf(stderr, "\t-v\t(v)erbose,       default is disabled, presents some debug information at runtime.\n\n");
}
