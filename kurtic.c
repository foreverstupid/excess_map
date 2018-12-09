#include "kurtic.h"

/*
 * Kurtic kernel
 */
inline static double kernel(double x, double s0, double s1)
{
    double xx = x * x;
    return exp(-0.5 * (s0 * xx + s1 * xx * xx) / (1 + xx));
}



/*
 * Gets excess and dispersion of kurtic kernel
 */
static void get_current_values(double *k, double *d, double s0, double s1,
    const struct vector_func *buf)
{
    double norm;
    double x;
    int i;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] = kernel(x, s0, s1);
        x += buf->grid.step;
    }

    norm = get_norm(buf);

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    *d = get_integral(buf) / norm;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    *k = get_integral(buf) / norm / (*d) / (*d) - 3;
}



int kurtic(const gsl_vector *x, void *params, gsl_vector *f)
{
    struct kurtic_params *p = (struct kurtic_params *)params;
    double curr_k;
    double curr_d;
    double s0 = gsl_vector_get(x, 0);
    double s1 = gsl_vector_get(x, 1);

    get_current_values(&curr_k, &curr_d, s0, s1, &(p->buffer));

    gsl_vector_set(f, 0, curr_k - p->k);
    gsl_vector_set(f, 1, curr_d - p->d);

    return GSL_SUCCESS;
}
