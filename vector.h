#ifndef VECTOR_MODULE_H
#define VECTOR_MODULE_H

#include <math.h>

/*
 * Linspace info
 */
struct linspace{
    int count;              /* step count */
    double step;            /* step value */
    double origin;          /* origin of space */
};




/*
 * A vector function representation
 */
struct vector_func{
    double *storage;                /* function values */
    struct linspace grid;           /* function grid */
};



/*
 * Calculates integral of given function
 */
double get_integral(const struct vector_func *f);

/*
 * Calculates an integral norm of given function
 */
double get_norm(const struct vector_func *f);

#endif
