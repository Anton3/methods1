#include <stdio.h>
#include <math.h>
#include "../methods/runge_kutta_system.h"

real f1_test1(real x, real u, real v) {
    return v - cosl(x);
}

real f2_test1(real x, real u, real v) {
    return u + sinl(x);
}

real u_test1(real x) {
    return -sinl(x);
}

real v_test1(real x) {
    return 0;
}

real f1_test2(real x, real u, real v) {
    return (2.1L)*v - u*u;
}

real f2_test2(real x, real u, real v) {
    return expl(-u) + x + (2.1L)*v;
}

real f1_test3(real x, real u, real v) {
    return (u - v) / x;
}

real f2_test3(real x, real u, real v) {
    return (u + v) / x;
}


const size_t BUF_SIZE = 80;

int main(void) {
    cauchy_problem_system problems[] = {
        // f1        f2        a  b   y0 z0    n
        {  f1_test1, f2_test1, 0, 10, 0, 0,    100  },
        {  f1_test1, f2_test1, 0, 10, 0, 0,    1000 },
        {  f1_test1, f2_test1, 0, 10, 0, 0,    3000 },
        {  f1_test2, f2_test2, 0, 1.5,1, 0.25, 300  },
        {  f1_test3, f2_test3, 1, 10, 1, 1,    100  },
    };

    char fname[BUF_SIZE];

    for (size_t i = 0; i < 5; ++i) {
        cauchy_solution_system solution = runge_kutta2_solve_system(problems[i]);

        snprintf(fname, BUF_SIZE, "runge-kutta2-system-%zu", i+1);
        print_cauchy_solution_system(solution, fname);

        if (i < 3) {
            printf("Runge-Kutta 2: test %zu, error u: %11.8Lf, error v: %11.8Lf\n", i+1,
                   cauchy_solution_error_u(solution, u_test1),
                   cauchy_solution_error_v(solution, v_test1));
        }
    }

    for (size_t i = 0; i < 5; ++i) {
        cauchy_solution_system solution = runge_kutta4_solve_system(problems[i]);

        snprintf(fname, BUF_SIZE, "runge-kutta4-system-%zu", i+1);
        print_cauchy_solution_system(solution, fname);

        if (i < 3) {
            printf("Runge-Kutta 4: test %zu, error u: %11.8Lf, error v: %11.8Lf\n", i+1,
                   cauchy_solution_error_u(solution, u_test1),
                   cauchy_solution_error_v(solution, v_test1));
        }
    }

    return 0;
}
