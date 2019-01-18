#include "solver.h"

/*
 * Initializes multidimensional solver.
 */
static gsl_multiroot_fdfsolver *make_solver()
{
    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_gnewton;
    return gsl_multiroot_fdfsolver_alloc(T, 2);
}



/*
 * Initializes result information
 */
static void init_result_info(struct result *res, struct problem_info p)
{
    int length = p.k_grid.count * p.d_grid.count;

    res->s0.storage = malloc(sizeof(double) * length);
    res->s0.grid.count = length;

    res->s1.storage = malloc(sizeof(double) * length);
    res->s1.grid.count = length;
}



/*
 * Solves an equation system
 */
static void find_root(
    double *s0,
    double *s1,
    gsl_multiroot_fdfsolver *solver,
    gsl_multiroot_function_fdf *f,
    int max_iter_count,
    double eps,
    double beg_s0,
    double beg_s1
)
{
    int status;
    size_t iter = 0;
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, beg_s0);
    gsl_vector_set(x, 1, beg_s1);

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

    *s0 = gsl_vector_get(solver->x, 0);
    *s1 = gsl_vector_get(solver->x, 1);

    printf("Solution: (s0 = %lf, s1 = %lf)\n\n", *s0, *s1);

    gsl_vector_free(x);
}



/*
 * Creates begin vector.
 */
void get_begin(double d, double k, double *s0, double *s1)
{
    if(d < 1 && k <= 0){
        *s0 = 0.1;
        *s1 = 0.02;
        return;
    }

    if(k > 0){
        *s0 = 10;
        *s1 = 1;
        return;
    }

    *s0 = 400;
    *s1 = 0.01;

#   ifdef DEBUG
    printf("Begin: (s0_0 = %lf, s1_0 = %lf)\n", *s0, *s1);
#   endif
}



/*
 * Solves the problem
 */
struct result solve(struct problem_info p)
{
    int i;
    int j;
    double k = p.k_grid.origin;
    double d;
    double s0;
    double s1;
    gsl_multiroot_fdfsolver *solver = make_solver();
    gsl_multiroot_function_fdf f = { &kurtic_f, &kurtic_df, &kurtic_fdf, 2, NULL };
    struct kurtic_params params;
    struct result res;

    init_result_info(&res, p);
    params.buffer.storage = malloc(sizeof(double) * p.space_grid.count);
    params.buffer.grid = p.space_grid;
    f.params = &params;

    for(i = 0; i < p.k_grid.count; i++){
        d = p.d_grid.origin;
        for(j = 0; j < p.d_grid.count; j++){
            params.k = k;
            params.d = d * d;

#           ifdef DEBUG
            printf("Space: [%lf; %lf]\n", params.buffer.grid.origin,
                params.buffer.grid.origin + params.buffer.grid.step *
                params.buffer.grid.count);
#           endif

            get_begin(params.d, params.k, &s0, &s1);
            printf("Input: (k = %lf, d = %lf)\n", k, params.d);

            find_root(
                res.s0.storage + i * p.d_grid.count + j,
                res.s1.storage + i * p.d_grid.count + j,
                solver,
                &f,
                p.iter_count,
                p.eps,
                s0,
                s1
            );

            d += p.d_grid.step;
        }

        k += p.k_grid.step;
    }

    free(params.buffer.storage);
    gsl_multiroot_fdfsolver_free(solver);

    return res;
}
