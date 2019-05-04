#include "solver.h"

/*
 * Initializes multidimensional fdf solver
 */
static gsl_multiroot_fdfsolver *make_fdf_solver()
{
    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_gnewton;
    return gsl_multiroot_fdfsolver_alloc(T, 2);
}



/*
 * Initializes multidimensional simple solver
 */
static gsl_multiroot_fsolver *make_f_solver()
{
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_dnewton;
    return gsl_multiroot_fsolver_alloc(T, 2);
}



/*
 * Initializes result information
 */
static void init_result_info(struct result *res, struct problem_info *p)
{
    int length = p->k_grid.count * p->d_grid.count;

    res->a.storage = malloc(sizeof(double) * length);
    res->a.grid.count = length;

    res->b.storage = malloc(sizeof(double) * length);
    res->b.grid.count = length;
}



/*
 * Creates begin solution vector
 */
static void get_begin(struct problem_info *p, double d, double k, double *a,
    double *b)
{
    if(p->kern_type == KURTIC){
        if(d < 1 && k <= 0){
            *a = 0.1;
            *b = 0.02;
            return;
        }

        if(k > 0){
            *a = 10;
            *b = 1;
            return;
        }

        *a = 400;
        *b = 0.01;
    }else if(p->kern_type == RGARDEN){
        *a = 1;
        if(k < 0){
            *b = 4;
        }else if(fabs(k) < 10e-7){
            *b = 2;
        }else{
            *b = 1;
        }
    }else if(p->kern_type == POLYEXP){
        if(k < 0){
            *b = 1;
            *a = 0;
        }else{
            *b = 0;
            *a = 1;
        }
    }else{
        *a = 1;
        *b = 1;
    }

#   ifdef DEBUG
    printf("Begin: (a_0 = %lf, b_0 = %lf)\n", *a, *b);
#   endif
}



/*
 * Solves an equation system using fdf solver
 */
static void find_root_fdf(
    double *a,
    double *b,
    gsl_multiroot_fdfsolver *solver,
    gsl_multiroot_function_fdf *f,
    int max_iter_count,
    double eps,
    double beg_a,
    double beg_b
)
{
    int status;
    size_t iter = 0;
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, beg_a);
    gsl_vector_set(x, 1, beg_b);

    gsl_multiroot_fdfsolver_set(solver, f, x);

    do{
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(solver);

        if(status){
            printf("Stucked! (%ld)\n", iter);
            break;
        }

        status = gsl_multiroot_test_residual(solver->f, eps);
    }while(status == GSL_CONTINUE && iter < max_iter_count);

    *a = gsl_vector_get(solver->x, 0);
    *b = gsl_vector_get(solver->x, 1);

    gsl_vector_free(x);
}



/*
 * Solves the problem using derevative method
 */
static struct result solve_fdf(struct problem_info *p)
{
    int i;
    int j;
    double k = p->k_grid.origin;
    double d;
    double a;
    double b;
    gsl_multiroot_fdfsolver *solver = make_fdf_solver();
    gsl_multiroot_function_fdf f = { p->f, p->df, p->fdf, 2, NULL };
    struct params params;
    struct result res;

    init_result_info(&res, p);
    params.buffer.storage = malloc(sizeof(double) * p->space_grid.count);
    params.buffer.grid = p->space_grid;
    f.params = &params;

    for(i = 0; i < p->k_grid.count; i++){
        d = p->d_grid.origin;
        for(j = 0; j < p->d_grid.count; j++){
            params.k = k;
            params.d = d * d;

#           ifdef DEBUG
            printf("Space: [%lf; %lf]\n", params.buffer.grid.origin,
                params.buffer.grid.origin + params.buffer.grid.step *
                params.buffer.grid.count);
#           endif

            get_begin(p, params.d, params.k, &a, &b); 
            find_root_fdf(
                res.a.storage + i * p->d_grid.count + j,
                res.b.storage + i * p->d_grid.count + j,
                solver,
                &f,
                p->iter_count,
                p->eps,
                a,
                b
            );
            printf(
                "Input: (k = %lf, d = %lf)\n"
                "Solution: (a= %lf, b= %lf)\n\n",
                k,
                params.d,
                res.a.storage[i * p->d_grid.count + j],
                res.b.storage[i * p->d_grid.count + j]
            );

            d += p->d_grid.step;
        }

        k += p->k_grid.step;
    }

    free(params.buffer.storage);
    gsl_multiroot_fdfsolver_free(solver);

    return res;
}



/*
 * Solves an equation system using simple iterative solver
 */
static void find_root_f(
    double *a,
    double *b,
    gsl_multiroot_fsolver *solver,
    gsl_multiroot_function *f,
    int max_iter_count,
    double eps,
    double beg_a,
    double beg_b
)
{
    int status;
    size_t iter = 0;
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, beg_a);
    gsl_vector_set(x, 1, beg_b);

    gsl_multiroot_fsolver_set(solver, f, x);

    do{
        iter++;
        status = gsl_multiroot_fsolver_iterate(solver);

        if(status){
            printf("Stucked! (%ld)\n", iter);
            break;
        }

        status = gsl_multiroot_test_residual(solver->f, eps);
    }while(status == GSL_CONTINUE && iter < max_iter_count);

    *a = gsl_vector_get(solver->x, 0);
    *b = gsl_vector_get(solver->x, 1);

    gsl_vector_free(x);
}



/*
 * Solves the problem using simple iterative method
 */
static struct result solve_f(struct problem_info *p)
{
    int i;
    int j;
    double k = p->k_grid.origin;
    double d;
    double a;
    double b;
    gsl_multiroot_fsolver *solver = make_f_solver();
    gsl_multiroot_function f = { p->f, 2, NULL };
    struct params params;
    struct result res;

    init_result_info(&res, p);
    params.buffer.storage = malloc(sizeof(double) * p->space_grid.count);
    params.buffer.grid = p->space_grid;
    f.params = &params;

    for(i = 0; i < p->k_grid.count; i++){
        d = p->d_grid.origin;
        for(j = 0; j < p->d_grid.count; j++){
            params.k = k;
            params.d = d * d;

#           ifdef DEBUG
            printf("Space: [%lf; %lf]\n", params.buffer.grid.origin,
                params.buffer.grid.origin + params.buffer.grid.step *
                params.buffer.grid.count);
#           endif

            get_begin(p, params.d, params.k, &a, &b); 
            find_root_f(
                res.a.storage + i * p->d_grid.count + j,
                res.b.storage + i * p->d_grid.count + j,
                solver,
                &f,
                p->iter_count,
                p->eps,
                a,
                b
            );
            printf(
                "Input: (k = %lf, d = %lf)\n"
                "Solution: (a = %lf, b = %lf)\n\n",
                k,
                params.d,
                res.a.storage[i * p->d_grid.count + j],
                res.b.storage[i * p->d_grid.count + j]
            );

            d += p->d_grid.step;
        }

        k += p->k_grid.step;
    }

    free(params.buffer.storage);
    gsl_multiroot_fsolver_free(solver);

    return res;
}



struct result solve(struct problem_info *p)
{
    if(p->kern_type == KURTIC){
        return solve_fdf(p);
    }

    return solve_f(p);
}
