#include <stdio.h>
#include "runge_kutta_system.h"

void print_cauchy_solution_system(cauchy_solution_system s, const char *fname) {
    FILE *output = fopen(fname, "w");
    if (output == NULL) { exit(1); }

    size_t n = s.x->dimension;

    for (size_t i = 0; i < n; ++i) {
        fprintf(output, "%-11.8Lf %-11.8Lf\n", vidx(s.x, i), vidx(s.y, i));
    }

    fprintf(output, "\n\n");

    for (size_t i = 0; i < n; ++i) {
        fprintf(output, "%-11.8Lf %-11.8Lf\n", vidx(s.x, i), vidx(s.z, i));
    }

    fclose(output);
}

real cauchy_solution_error_u(cauchy_solution_system s, one_arg_func u) {
    size_t n = s.x->dimension;
    vector actual_y = new_vector(n);

    for (size_t i = 0; i < n; ++i) {
        vidx(actual_y, i) = u(vidx(s.x, i));
    }

    real result = vector_distance(s.y, actual_y);
    delete_vector(actual_y);
    return result;
}

real cauchy_solution_error_v(cauchy_solution_system s, one_arg_func v) {
    size_t n = s.x->dimension;
    vector actual_z = new_vector(n);

    for (size_t i = 0; i < n; ++i) {
        vidx(actual_z, i) = v(vidx(s.x, i));
    }

    real result = vector_distance(s.z, actual_z);
    delete_vector(actual_z);
    return result;
}

cauchy_solution_system runge_kutta2_solve_system(cauchy_problem_system p) {
    real h = (p.b - p.a) / p.n;

    vector x = new_vector(p.n + 1);
    vidx(x, 0) = p.a;

    vector y = new_vector(p.n + 1);
    vidx(y, 0) = p.y0;

    vector z = new_vector(p.n + 1);
    vidx(z, 0) = p.z0;

    real x_i = p.a;
    real y_i = p.y0;
    real z_i = p.z0;

    for (size_t i = 0; i < p.n; ++i) {
        real x_ip1 = ((p.n-(i+1)) * p.a + (i+1) * p.b) / p.n;

        real k1 = p.f1(x_i, y_i, z_i);
        real l1 = p.f2(x_i, y_i, z_i);

        real k2 = p.f1(x_ip1, y_i + h*k1, z_i + h*l1);
        real l2 = p.f2(x_ip1, y_i + h*k1, z_i + h*l1);

        real y_ip1 = y_i + (h/2) * (k1 + k2);
        real z_ip1 = z_i + (h/2) * (l1 + l2);

        vidx(x, i+1) = x_i = x_ip1;
        vidx(y, i+1) = y_i = y_ip1;
        vidx(z, i+1) = z_i = z_ip1;
    }

    cauchy_solution_system solution = { x, y, z };
    return solution;
}

cauchy_solution_system runge_kutta4_solve_system(cauchy_problem_system p) {
    real h = (p.b - p.a) / p.n;

    vector x = new_vector(p.n + 1);
    vidx(x, 0) = p.a;

    vector y = new_vector(p.n + 1);
    vidx(y, 0) = p.y0;

    vector z = new_vector(p.n + 1);
    vidx(z, 0) = p.z0;

    real x_i = p.a;
    real y_i = p.y0;
    real z_i = p.z0;

    for (size_t i = 0; i < p.n; ++i) {
        real x_ip1 = ((p.n-(i+1)) * p.a + (i+1) * p.b) / p.n;
        real x_half = x_i + h/2;

        real k1 = p.f1(x_i, y_i, z_i);
        real l1 = p.f2(x_i, y_i, z_i);

        real k2 = p.f1(x_half, y_i + (h/2)*k1, z_i + (h/2)*l1);
        real l2 = p.f2(x_half, y_i + (h/2)*k1, z_i + (h/2)*l1);

        real k3 = p.f1(x_half, y_i + (h/2)*k2, z_i + (h/2)*l2);
        real l3 = p.f2(x_half, y_i + (h/2)*k2, z_i + (h/2)*l2);

        real k4 = p.f1(x_ip1, y_i + h*k3, z_i + h*l3);
        real l4 = p.f2(x_ip1, y_i + h*k3, z_i + h*l3);

        real y_ip1 = y_i + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        real z_ip1 = z_i + (h/6) * (l1 + 2*l2 + 2*l3 + l4);

        vidx(x, i+1) = x_i = x_ip1;
        vidx(y, i+1) = y_i = y_ip1;
        vidx(z, i+1) = z_i = z_ip1;
    }

    cauchy_solution_system solution = { x, y, z };
    return solution;
}
