#include <stdio.h>
#include <math.h>
#include "../base/vectors.h"
#include "../methods/runge_kutta.h"

real f_test1(real x, real y) {
    return 3 - y - x;
}

real u_test1(real x) {
    return 4 - x - 4*expl(-x);
}

real f_test2(real x, real y) {
    return sinl(x) - y;
}

real u_test2(real x) {
    return 0.5 * (sin(x) - cos(x) + 21*expl(-x));
}

real f_test3(real x, real y) {
    return -y - x*x;
}

real u_test3(real x) {
    return (2-x)*x - 2 + 12*expl(-x);
}


const size_t BUF_SIZE = 80;

int main(void) {
    cauchy_problem problems[] = {
        // f        a  b   y0 n
        {  f_test1, 0, 5,  0,  10   },
        {  f_test1, 0, 5,  0,  100  },
        {  f_test1, 0, 5,  0,  1000 },
        {  f_test2, 0, 10, 10, 100  },
        {  f_test3, 0, 10, 10, 100  }
    };

    one_arg_func checks[] = {
        u_test1,
        u_test1,
        u_test1,
        u_test2,
        u_test3
    };

    char fname[BUF_SIZE];

    for (size_t i = 0; i < 5; ++i) {
        cauchy_solution solution = runge_kutta2_solve(problems[i]);

        snprintf(fname, BUF_SIZE, "runge-kutta2-%zu", i+1);
        print_cauchy_solution(solution, fname);

        printf("Runge-Kutta 2: test %zu, error: %11.8Lf\n", i+1,
               cauchy_solution_error(solution, checks[i]));
    }

    for (size_t i = 0; i < 5; ++i) {
        cauchy_solution solution = runge_kutta4_solve(problems[i]);

        snprintf(fname, BUF_SIZE, "runge-kutta4-%zu", i+1);
        print_cauchy_solution(solution, fname);

        printf("Runge-Kutta 4: test %zu, error: %11.8Lf\n", i+1,
               cauchy_solution_error(solution, checks[i]));
    }

    return 0;
}
