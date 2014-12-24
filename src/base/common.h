#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>

typedef long double real;

typedef real (*one_arg_func)(real);
typedef real (*two_arg_func)(real, real);
typedef real (*three_arg_func)(real, real, real);

#define swap(a, b) do { typeof(a) temp = a;                            \
                        a = b;                                         \
                        b = temp;                                      \
                   } while (0)

bool is_zero(real value);

extern const int OUTPUT_WIDTH;

#endif
