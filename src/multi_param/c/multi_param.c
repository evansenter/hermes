#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "multi_param.h"

PARAM_CONTAINER* split_args(int argc, char** argv, char** subparam_order, int subparam_count) {
  int i, last_subarg = -1;
  PARAM_CONTAINER* param_container;

  SUBARG_MATCH current_match = {
    .argv  = NULL,
    .index = -1
  };

  param_container = malloc(subparam_count * sizeof(PARAM_CONTAINER));
  for (i = 0; i < subparam_count; ++i) {
    param_container[i].argv = malloc(argc * sizeof(char*));
    param_container[i].argc = 0;
  }

  for (i = 1; i < argc; ++i) {
    match_for_subparam(subparam_order, argv[i], subparam_count, &current_match);
    if (current_match.index != -1) {
      last_subarg = current_match.index;
      param_container[last_subarg].argv[param_container[last_subarg].argc++] = current_match.argv;

      #ifdef DEBUG
        printf("Updated last_subarg to %d\n", last_subarg);
      #endif
    } else if (last_subarg >= 0) {
      param_container[last_subarg].argv[param_container[last_subarg].argc++] = argv[i];

      #ifdef DEBUG
        printf("Using last_subarg position %d\n", last_subarg);
      #endif
    }
  }

  return param_container;
}

void match_for_subparam(char** subparam_order, char* subparam, int subparam_count, SUBARG_MATCH* match) {
  int i;
  char* token;
  token = strtok(subparam, "-");

  while (token != NULL) {
    for (i = 0; i < subparam_count; ++i) {
      if (!strcmp(subparam_order[i], token)) {
        token = strtok(NULL, "-");
        match->argv  = strdup(token);
        match->index = i;

        #ifdef DEBUG
          printf("Found location for %s at index %d\n", match->argv, match->index);
        #endif

        return;
      }
    }

    token = strtok(NULL, "-");
  }

  #ifdef DEBUG
    printf("Couldn't find location for %s\n", subparam);
  #endif

  match->index = -1;
}

void print_multi_params(PARAM_CONTAINER* params, char** subparam_order, int subparam_count) {
  int i, j;

  for (i = 0; i < subparam_count; ++i) {
    if (params[i].argc > 0) {
      printf("Params for %s:\n\t", subparam_order[i]);

      for (j = 0; j < params[i].argc; ++j) {
        printf("%s ", params[i].argv[j]);
      }

      printf("\n");
    }
  }
}
