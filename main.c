#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "vector.h"
#include "solver.h"

/*
 * Holds info about writing output data
 */
struct output_info{
    const char *file_name;      /* output file name */
    struct linspace k_grid;     /* excess grid from problem info */
    struct linspace d_grid;     /* dispersion grid from problem info */
};



/*
 * Initializes problem info
 */
void make_problem_info(int argc, const char **argv, struct problem_info **p)
{
    double beg;
    int count;
    double last;

    sscanf(argv[1], "%lf", &(beg));
    sscanf(argv[2], "%d", &(count));
    sscanf(argv[3], "%lf", &(last));
    (*p)->k_grid.origin = beg;
    (*p)->k_grid.count = count;
    (*p)->k_grid.step = count > 1 ? (last - beg) / (count - 1) : 0.0;

    sscanf(argv[4], "%lf", &(beg));
    sscanf(argv[5], "%d", &(count));
    sscanf(argv[6], "%lf", &(last));
    (*p)->d_grid.origin = beg;
    (*p)->d_grid.count = count;
    (*p)->d_grid.step = count > 1 ? (last - beg) / (count - 1) : 0.0;

    sscanf(argv[7], "%d", &(count));
    (*p)->space_grid.count = count;

    sscanf(argv[8], "%lf", &((*p)->eps));
    (*p)->iter_count = 100;

    (*p)->kern_type = argv[9][0];
    if((*p)->kern_type == KURTIC){
        (*p)->f = kurtic_f;
        (*p)->df = kurtic_df;
        (*p)->fdf = kurtic_fdf;
    }else if((*p)->kern_type == RGARDEN){
        (*p)->f = rgarden_f;
    }else{
        *p = NULL;
    }
}



/*
 * Initializes output data writing info
 */
struct output_info make_output_info(int argc, const char **argv,
    const struct problem_info *p)
{
    struct output_info info = { argv[10], p->k_grid, p->d_grid };
    return info;
}



/*
 * Writes solution
 */
void print(struct result res, struct output_info oinf)
{
    FILE *out = fopen(oinf.file_name, "w");
    int i;
    int j;
    int index;
    double k = oinf.k_grid.origin;
    double d = oinf.d_grid.origin;

    for(i = 0; i < oinf.k_grid.count; i++){
        index = i * oinf.d_grid.count;
        d = oinf.d_grid.origin;
        for(j = 0; j < oinf.d_grid.count; j++){
            fprintf(
                out,
                "%lf %lf %lf %lf\n",
                k,
                d,
                res.a.storage[index + j],
                res.b.storage[index + j]
            );

            d += oinf.d_grid.step;
        }

        k += oinf.k_grid.step;
    }

    fclose(out);
}



#ifdef DEBUG
/*
 * Prints information get from cmd
 */
void print_given_info(struct problem_info *p, struct output_info o)
{
    printf("Excess kurtosis grid:\n");
    printf("    origin: %lf\n", p->k_grid.origin);
    printf("      step: %lf\n", p->k_grid.step);
    printf("     count: %d\n", p->k_grid.count);

    printf("\nDispersion grid:\n");
    printf("    origin: %lf\n", p->d_grid.origin);
    printf("      step: %lf\n", p->d_grid.step);
    printf("     count: %d\n", p->d_grid.count);

    printf("\nSpace grid:\n");
    printf("    origin: %lf\n", p->space_grid.origin);
    printf("      step: %lf\n", p->space_grid.step);
    printf("     count: %d\n", p->space_grid.count);

    printf("\nOutput file: %s\n", o.file_name);
}
#endif



int main(int argc, const char **argv)
{
    struct problem_info *prinf = malloc(sizeof(struct problem_info));
    struct output_info oinf;
    struct result res;
    
    make_problem_info(argc, argv, &prinf);
    if(prinf == NULL){
        fprintf(stderr, "### Invalid arguments!\n");
        return 1;
    }

#   ifdef DEBUG
    print_given_info(prinf, oinf);
#   endif
    
    oinf = make_output_info(argc, argv, prinf);
    res = solve(prinf);
    print(res, oinf);

    free(res.a.storage);
    free(res.b.storage);
    free(prinf);

    return 0;
}
