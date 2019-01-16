#include "kurtic.h"

/*
 * Kurtic kernel
 */
inline static double kernel(double x, double s0, double s1)
{
    double xx = x * x;
    return exp(-0.5 * (s0 * xx + s1 * xx * xx) / (1 + xx));
}



static double get_origin(double s0, double s1)
{
    double x = 0.0;
    double step = 1e-5;

    while(fabs(kernel(x, s0, s1)) > 1e-12){
        x -= step;
    }

#   ifdef DEBUG
    printf("Origin: %lf\n", x);
#   endif

    return x;
}



static void get_current_values_df(double s0, double s1, struct vector_func *buf, gsl_matrix *J)
{
    double norm;
    double mu;
    double sgm;
    double dsgms0;
    double dsgms1;
    double dmus0;
    double dmus1;
    double x;
    double xx;
    double d;
    double dd;
    double d4;
    double tmp;
    int i;

    buf->grid.origin = get_origin(s0, s1);
    buf->grid.step = 2 * fabs(buf->grid.origin) / (buf->grid.count - 1);

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

    sgm = get_integral(buf);
    d = sgm / norm;
    dd = d * d;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    mu = get_integral(buf) / norm;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        xx = x * x;
        buf->storage[i] /= 1 + xx;
        x += buf->grid.step;
    }

    dsgms0 = -0.5 * get_integral(buf) / norm;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    dsgms1 = -0.5 * get_integral(buf) / norm;

    d4 = dd * dd;
    dmus0 = dsgms1 / dd - 2 * mu * d * dsgms0 / d4;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    tmp = -0.5 * get_integral(buf) / norm;
    dmus1 = tmp / dd - 2 * mu * d * dsgms1 / d4;

    gsl_matrix_set(J, 0, 0, dmus0);
    gsl_matrix_set(J, 0, 1, dmus1);
    gsl_matrix_set(J, 1, 0, dsgms0);
    gsl_matrix_set(J, 1, 1, dsgms1);

#   ifdef DEBUG
    printf("J = [ %10.3lf, %10.3lf ]\n", dmus0, dmus1);
    printf("    [ %10.3lf, %10.3lf ]\n", dsgms0, dsgms1);
#   endif
}





static void get_current_values_fdf(double *k, double *d, double s0, double s1,
    struct vector_func *buf, gsl_matrix *J)
{
    double norm;
    double mu;
    double sgm;
    double dsgms0;
    double dsgms1;
    double dmus0;
    double dmus1;
    double x;
    double xx;
    double dd;
    double d4;
    double tmp;
    int i;

    buf->grid.origin = get_origin(s0, s1);
    buf->grid.step = 2 * fabs(buf->grid.origin) / (buf->grid.count - 1);

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

    sgm = get_integral(buf);
    *d = sgm / norm;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    mu = get_integral(buf) / norm;
    dd = (*d) * (*d);
    *k = mu / dd - 3;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        xx = x * x;
        buf->storage[i] /= 1 + xx;
        x += buf->grid.step;
    }

    dsgms0 = -0.5 * get_integral(buf) / norm;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    dsgms1 = -0.5 * get_integral(buf) / norm;

    d4 = dd * dd;
    dmus0 = dsgms1 / dd - 2 * mu * (*d) * dsgms0 / d4;

    x = buf->grid.origin;
    for(i = 0; i < buf->grid.count; i++){
        buf->storage[i] *= x * x;
        x += buf->grid.step;
    }

    tmp = -0.5 * get_integral(buf) / norm;
    dmus1 = tmp / dd - 2 * mu * (*d) * dsgms1 / d4;

    gsl_matrix_set(J, 0, 0, dmus0);
    gsl_matrix_set(J, 0, 1, dmus1);
    gsl_matrix_set(J, 1, 0, dsgms0);
    gsl_matrix_set(J, 1, 1, dsgms1);
}



static void get_current_values_f(double *k, double *d, double s0, double s1,
    struct vector_func *buf)
{
    double norm;
    double x;
    int i;

    buf->grid.origin = get_origin(s0, s1);
    buf->grid.step = 2 * fabs(buf->grid.origin) / (buf->grid.count - 1);

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
#       ifdef DEBUG
        printf("k = %lf, d = %lf\n", *k, *d);
        printf("{s0 = %lf, s1 = %lf}\n", s0, s1);
#       endif
}

inline static double psgn(double x)
{
    return 0.0;//x > 0 ? 0.0 : 1.0;
}



int kurtic_fdf(const gsl_vector *x, void *params, gsl_vector *f,
    gsl_matrix *J)
{
    struct kurtic_params *p = (struct kurtic_params *)params;
    double curr_k;
    double curr_d;
    double s0 = gsl_vector_get(x, 0);
    double s1 = gsl_vector_get(x, 1);

    get_current_values_fdf(&curr_k, &curr_d, s0, s1, &(p->buffer), J);

    gsl_vector_set(f, 0, curr_k - p->k + psgn(s0) + psgn(s1));
    gsl_vector_set(f, 1, curr_d - p->d + psgn(s0) + psgn(s1));

    return GSL_SUCCESS;
}



int kurtic_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    struct kurtic_params *p = (struct kurtic_params *)params;
    double s0 = gsl_vector_get(x, 0);
    double s1 = gsl_vector_get(x, 1);

    get_current_values_df(s0, s1, &(p->buffer), J);

    return GSL_SUCCESS;
}



int kurtic_f(const gsl_vector *x, void *params, gsl_vector *f)
{
    struct kurtic_params *p = (struct kurtic_params *)params;
    double curr_k;
    double curr_d;
    double s0 = gsl_vector_get(x, 0);
    double s1 = gsl_vector_get(x, 1);

    get_current_values_f(&curr_k, &curr_d, s0, s1, &(p->buffer));

    gsl_vector_set(f, 0, curr_k - p->k + psgn(s0) + psgn(s1));
    gsl_vector_set(f, 1, curr_d - p->d + psgn(s0) + psgn(s1));

    return GSL_SUCCESS;
}
