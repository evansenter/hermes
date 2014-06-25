#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "multi_param.h"

PARAM_CONTAINER* split_args(int argc, char** argv, char** subparam_order, int subparam_count) {
  int i, j, last_subarg = NO_MATCH;
  PARAM_CONTAINER* param_container;

  SUBARG_MATCH current_match = {
    .argv  = NULL,
    .index = -1
  };

  param_container = malloc(subparam_count * sizeof(PARAM_CONTAINER));

  for (i = 0; i < subparam_count; ++i) {
    param_container[i].argv    = malloc(argc * sizeof(char*));
    param_container[i].argv[0] = strdup(subparam_order[i]);
    param_container[i].argc    = 1;
  }

  for (i = 1; i < argc; ++i) {
    match_for_subparam(subparam_order, argv[i], subparam_count, &current_match);

    if (current_match.index == ALL_SUBARG) {
      // If the current match is prefixed --all or the current match isn't prefixed and the last flag was --all (meaning it's a value).
      last_subarg = ALL_SUBARG;

      for (j = 0; j < subparam_count; ++j) {
        param_container[j].argv[param_container[j].argc++] = current_match.argv;
      }

#ifdef DEBUG
      printf("Updated last_subarg to %d\n", last_subarg);
#endif
    } else if (current_match.index == NO_MATCH && last_subarg == ALL_SUBARG) {
      // It's not a global flag and it isn't prefixed (meaning it's a value for a global flag).

      for (j = 0; j < subparam_count; ++j) {
        param_container[j].argv[param_container[j].argc++] = argv[i];
      }
    } else if (current_match.index != NO_MATCH) {
      // It's not a global flag and we can identify the prefix.
      last_subarg = current_match.index;
      param_container[last_subarg].argv[param_container[last_subarg].argc++] = current_match.argv;

#ifdef DEBUG
      printf("Updated last_subarg to %d\n", last_subarg);
#endif
    } else if (last_subarg >= 0) {
      // It's not a global flag and it isn't prefixed (meaning it's a value for a specific subparam).
      param_container[last_subarg].argv[param_container[last_subarg].argc++] = argv[i];
    }
  }

  return param_container;
}

void match_for_subparam(char** subparam_order, char* subparam, int subparam_count, SUBARG_MATCH* match) {
  int i;
  char* token;
  char* token_with_prefix;

  if (!strncmp(subparam, "--", 2)) {
    token = strtok(subparam, "-");

    while (token != NULL) {
      if (!strcmp("all", token)) {
        token             = strtok(NULL, "-");
        token_with_prefix = malloc(strlen(token) + 2);
        strcpy(token_with_prefix, "-");
        strcat(token_with_prefix, token);
        match->argv       = token_with_prefix;
        match->index      = ALL_SUBARG;

#ifdef DEBUG
        printf("Found global flag for %s\n", match->argv);
#endif
        return;
      } else {
        for (i = 0; i < subparam_count; ++i) {
          if (!strcmp(subparam_order[i], token)) {
            token             = strtok(NULL, "-");
            token_with_prefix = malloc(strlen(token) + 2);
            strcpy(token_with_prefix, "-");
            strcat(token_with_prefix, token);
            match->argv       = token_with_prefix;
            match->index      = i;

#ifdef DEBUG
            printf("Found location for %s at index %d\n", match->argv, match->index);
#endif
            return;
          }
        }
      }

      token = strtok(NULL, "-");
    }
  }

#ifdef DEBUG
  printf("Couldn't find location for %s\n", subparam);
#endif

  match->index = NO_MATCH;
}

void print_available_subflags(char* program, char* prefix, void (*flags)()) {
  fprintf(stderr, "Available %s flags (prefixed with --%s-) are:\n\n", program, prefix);
  (*flags)();
}

void print_valid_multi_params_prefixes(char** subparam_order, int subparam_count) {
  int i;

  fprintf(stderr, "This program supports the following flag prefixes:\n");

  fprintf(stderr, "\t--all- (i.e. --all-v for all Hermes libraries to use verbose mode)\n");
  for (i = 0; i < subparam_count; ++i) {
      fprintf(stderr, "\t--%s-\n", subparam_order[i]);
  }

  fprintf(stderr, "\n");
}

void print_multi_params(PARAM_CONTAINER* params, char** subparam_order, int subparam_count) {
  int i, j;

  for (i = 0; i < subparam_count; ++i) {
    if (params[i].argc > 0) {
      printf("Params for %s:\n", subparam_order[i]);

      for (j = 0; j < params[i].argc; ++j) {
        printf("\t%s\n", params[i].argv[j]);
      }

      printf("\n");
    }
  }
}
