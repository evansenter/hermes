#ifndef SHARED_TIMERS_H
#define SHARED_TIMERS_H

#include <sys/time.h>

#define START_ALL_TIMING                             \
  struct timeval full_start, full_stop, start, stop; \
  if (parameters.benchmark) {                        \
    gettimeofday(&full_start, NULL);                 \
    gettimeofday(&start, NULL);                      \
  }
  
#define STOP_ALL_TIMING                    \
  if (parameters.benchmark) {              \
    gettimeofday(&full_stop, NULL);        \
    PRINT_TIMING(full_start, full_stop, "total") \
  }

#define TIMING(task)                 \
  if (parameters.benchmark) {     \
    gettimeofday(&stop, NULL);       \
    PRINT_TIMING(start, stop, task); \
    gettimeofday(&start, NULL);      \
  }
  
#define PRINT_TIMING(start, stop, task) \
  printf("[benchmarking] %8.2f\ttime in ms for %s\n", (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0), task);

#define STOP_TIMING(task)           \
  if (parameters.benchmark) {       \
    gettimeofday(&stop, NULL);      \
    PRINT_TIMING(start, stop, task) \
  }

#endif
