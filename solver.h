#ifndef SOLVER_MODULE_H
#define SOLVER_MODULE_H

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "kurtic.h"
#include "vector.h"

/*
 * Holds info about problem initial data
 */
struct problem_info{
    struct linspace k_grid;         /* grid of kurtosis excess */
    struct linspace d_grid;         /* grid of dispersion */

    struct linspace space_grid;     /* grid of space */
    int iter_count;
    double eps;
};



/*
 * Holds info about problem solution (kurtic parameters)
 */
struct result{
    struct vector_func s0;
    struct vector_func s1;
};



/*
 * Solves the problem
 */
struct result solve(struct problem_info p);

#endif
