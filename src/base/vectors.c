#include <math.h>
#include "vectors.h"

vector new_vector(size_t n) {
    vector result = malloc(sizeof(struct vector_m));
    result->storage = malloc(n * sizeof(real));
    result->dimension = n;
    return result;
}

void delete_vector(vector vec) {
    free(vec->storage);
    free(vec);
}

real vector_distance(vector v1, vector v2) {
    size_t n = v1->dimension;
    real max = 0;

    for (size_t i = 0; i < n; ++i) {
        real next = fabsl(vidx(v1, i) - vidx(v2, i));
        if (max < next) { max = next; }
    }

    return max;
}