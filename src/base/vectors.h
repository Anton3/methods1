#ifndef VECTORS_H
#define VECTORS_H

#include <stdlib.h>
#include "common.h"

typedef struct vector_m {
	real *storage;
	size_t dimension;
} *vector;

vector new_vector(size_t n);
void delete_vector(vector vec);

real vector_distance(vector v1, vector v2);

#define vidx(vec, idx) vec->storage[idx]

#endif
