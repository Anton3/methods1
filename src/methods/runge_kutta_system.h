#include "../base/common.h"
#include "../base/vectors.h"

typedef struct cauchy_problem_system {
    three_arg_func f1;
    three_arg_func f2;
    real a, b;
    real y0;
    real z0;
    size_t n;
} cauchy_problem_system;

typedef struct cauchy_solution_system {
    vector x;
    vector y;
    vector z;
} cauchy_solution_system;

void print_cauchy_solution_system(cauchy_solution_system solution, const char *fname);
real cauchy_solution_error_u(cauchy_solution_system s, one_arg_func u);
real cauchy_solution_error_v(cauchy_solution_system s, one_arg_func v);

cauchy_solution_system runge_kutta2_solve_system(cauchy_problem_system p);
cauchy_solution_system runge_kutta4_solve_system(cauchy_problem_system p);