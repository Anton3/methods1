#include <math.h>
#include "common.h"

const int OUTPUT_WIDTH = 18;
const real ZERO_EPSILON = 1e-50;

bool is_zero(real value) {
    return fabsl(value) < ZERO_EPSILON;
}