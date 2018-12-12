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
    double eps
)
{
    int status;
    size_t iter = 0;
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, 5.0);
    gsl_vector_set(x, 1, 1.0);

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

    printf("Solution: (%lf, %lf)\n", *s0, *s1);

    gsl_vector_free(x);
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
            params.d = d;

            find_root(
                res.s0.storage + i * p.k_grid.count + j,
                res.s1.storage + i * p.d_grid.count + j,
                solver,
                &f,
                p.iter_count,
                p.eps
            );

            d += p.d_grid.step;
        }

        k += p.k_grid.step;
    }

    free(params.buffer.storage);
    gsl_multiroot_fdfsolver_free(solver);

    return res;
}