#ifndef MULTI_PARAM_H
#define MULTI_PARAM_H

typedef struct {
  char** argv;
  int argc;
} PARAM_CONTAINER;

typedef struct {
  char* argv;
  int index;
} SUBARG_MATCH;

// #define DEBUG 1

PARAM_CONTAINER* split_args(int, char**, char**, int);
void match_for_subparam(char**, char*, int, SUBARG_MATCH*);
void print_multi_params(PARAM_CONTAINER*, char**, int);

#endif
