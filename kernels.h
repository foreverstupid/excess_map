#ifndef KERNELS_MODULE_H
#define KERNELS_MODULE_H

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "vector.h"

/*
 * Params for calculation method
 */
struct params{
    double k;                   /* excess kurtosis value */
    double d;                   /* dispersion value */

    struct vector_func buffer;  /* calculation buffer */
};


/* Kurtic kernel */
int kurtic_f(const gsl_vector *x, void *params, gsl_vector *f);

int kurtic_df(const gsl_vector *x, void *params, gsl_matrix *J);

int kurtic_fdf(const gsl_vector *x, void *params, gsl_vector *f,
    gsl_matrix *J);



/* Exonent polynomial kernel */
int polyexp_f(const gsl_vector *x, void *params, gsl_vector *f);

/*int kurtic_df(const gsl_vector *x, void *params, gsl_matrix *J);

int kurtic_fdf(const gsl_vector *x, void *params, gsl_vector *f,
    gsl_matrix *J);
*/


/* Roughgarden kernel */
int rgarden_f(const gsl_vector *x, void *params, gsl_vector *f);

#endif
