#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H
#include <stdio.h>

// ANSI/IEEE 754-1985 floating point precision
#define FLOATING_POINT_PRECISION                (1.e-16)

// Minimal value for division
#define FLOATING_POINT_MINIMAL_FOR_DIVISION     (1.e-64)
// Minimal value for comparison
#define FLOATING_POINT_MINIMAL_FOR_COMPARE      (1.e-12)

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define PRINT_LEN 10

#endif // COMMON_DEFS_H

