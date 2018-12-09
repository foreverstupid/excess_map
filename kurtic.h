#ifndef KURTIC_MODULE_HPP
#define KURTIC_MODULE_HPP

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "vector.h"

/*
 * Params for kurtic functions
 */
struct kurtic_params{
    double k;                   /* excess kurtosis value */
    double d;                   /* dispersion value */

    struct vector_func buffer;  /* calculation buffer */
};



/*
 * Kurtic function calculator.
 * Calculates two values: dispersion and excess
 */
int kurtic(const gsl_vector *x, void *params, gsl_vector *f);

#endif
