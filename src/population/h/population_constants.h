#ifndef CONSTANTS_H
#define CONSTANTS_H

// #define DEBUG 1
// #define SUPER_HEAVY_DEBUG 1
// #define INSANE_DEBUG 1

#define RT (1e-3 * 1.9872041 * (273.15 + (temperature)))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define POS_INF(parameters) ((parameters).end_time + 1)
#define NEG_INF(parameters) ((parameters).start_time - 1)

#endif
