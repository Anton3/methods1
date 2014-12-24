#include "../base/common.h"
#include "../base/vectors.h"

typedef struct cauchy_problem {
    two_arg_func f;
    real a, b;
    real y0;
    size_t n;
} cauchy_problem;

typedef struct cauchy_solution {
    vector x;
    vector y;
} cauchy_solution;

void print_cauchy_solution(cauchy_solution solution, const char *fname);
real cauchy_solution_error(cauchy_solution solution, one_arg_func u);

cauchy_solution runge_kutta2_solve(cauchy_problem p);
cauchy_solution runge_kutta4_solve(cauchy_problem p);