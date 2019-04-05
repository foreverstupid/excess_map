#include "solver.h"

/*
 * Initializes multidimensional solver.
 */
static gsl_multiroot_fdfsolver *make_fdf_solver()
{
    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_gnewton;
    return gsl_multiroot_fdfsolver_alloc(T, 2);
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
 * Solves an equation system
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
 * Creates begin vector
 */
void get_begin(struct problem_info *p, double d, double k, double *a,
    double *b)
{
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

#   ifdef DEBUG
    printf("Begin: (s0_0 = %lf, s1_0 = %lf)\n", *a, *b);
#   endif
}



/*
 * Solves the problem using derevative method
 */
struct result solve_fdf(struct problem_info *p)
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
                "Solution: (s0 = %lf, s1 = %lf)\n\n",
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
