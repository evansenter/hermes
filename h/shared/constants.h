#ifndef SHARED_CONSTANTS_H
#define SHARED_CONSTANTS_H

#define INPUT_DEBUG 1
// #define DEBUG 1
// #define SUPER_HEAVY_DEBUG 1
// #define INSANE_DEBUG 1

#define VIENNA_RT (1e-3 * 1.9872041 * (273.15 + (temperature)))
#define RT (1e-3 * 1.9872041 * (273.15 + 37))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#endif
