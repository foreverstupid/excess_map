#ifndef SOLVER_MODULE_H
#define SOLVER_MODULE_H

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "kernels.h"
#include "vector.h"

typedef int (*FFunc)(const gsl_vector *, void *, gsl_vector *);
typedef int (*DFunc)(const gsl_vector *, void *, gsl_matrix *);
typedef int (*FDFunc)(const gsl_vector *, void *, gsl_vector *,
    gsl_matrix *);

/*
 * Holds info about problem initial data
 */
struct problem_info{
    struct linspace k_grid;         /* grid of kurtosis excess */
    struct linspace d_grid;         /* grid of dispersion */

    struct linspace space_grid;     /* grid of space */
    int iter_count;                 /* iteration max count */
    double eps;                     /* precision */

    FFunc f;                        /* function GSL representation */
    DFunc df;                       /* derivative GSL representation */
    FDFunc fdf;                     /* function and derivative GSL */
};



/*
 * Holds info about problem solution (kernel parameters)
 */
struct result{
    struct vector_func a;
    struct vector_func b;
};



/*
 * Solves the problem using derivative method
 */
struct result solve_fdf(struct problem_info *p);

#endif
