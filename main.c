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
struct problem_info make_problem_info(int argc, const char **argv)
{
    double beg;
    int count;
    double last;
    struct problem_info info;

    sscanf(argv[1], "%lf", &(beg));
    sscanf(argv[2], "%d", &(count));
    sscanf(argv[3], "%lf", &(last));
    info.k_grid.origin = beg;
    info.k_grid.count = count;
    info.k_grid.step = (last - beg) / (count - 1);

    sscanf(argv[4], "%lf", &(beg));
    sscanf(argv[5], "%d", &(count));
    sscanf(argv[6], "%lf", &(last));
    info.d_grid.origin = beg;
    info.d_grid.count = count;
    info.d_grid.step = (last - beg) / (count - 1);

    sscanf(argv[7], "%lf", &(beg));
    sscanf(argv[8], "%d", &(count));
    sscanf(argv[9], "%lf", &(last));
    info.space_grid.origin = beg;
    info.space_grid.count = count;
    info.space_grid.step = (last - beg) / (count - 1);

    sscanf(argv[10], "%lf", &(info.eps));
    info.iter_count = 100;

    return info;
}



/*
 * Initializes output data writing info
 */
struct output_info make_output_info(int argc, const char **argv,
    struct problem_info p)
{
    struct output_info info = { argv[11], p.k_grid, p.d_grid };
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
        for(j = 0; j < oinf.d_grid.count; j++){
            fprintf(
                out,
                "%lf %lf %lf %lf\n",
                k,
                d,
                res.s0.storage[index + j],
                res.s1.storage[index + j]
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
void print_given_info(struct problem_info p, struct output_info o)
{
    printf("Excess kurtosis grid:\n");
    printf("    origin: %lf\n", p.k_grid.origin);
    printf("      step: %lf\n", p.k_grid.step);
    printf("     count: %d\n", p.k_grid.count);

    printf("\nDispersion grid:\n");
    printf("    origin: %lf\n", p.d_grid.origin);
    printf("      step: %lf\n", p.d_grid.step);
    printf("     count: %d\n", p.d_grid.count);

    printf("\nSpace grid:\n");
    printf("    origin: %lf\n", p.space_grid.origin);
    printf("      step: %lf\n", p.space_grid.step);
    printf("     count: %d\n", p.space_grid.count);

    printf("\nOutput file: %s\n", o.file_name);
}
#endif



int main(int argc, const char **argv)
{
    struct problem_info prinf = make_problem_info(argc, argv);
    struct output_info oinf = make_output_info(argc, argv, prinf);

#   ifdef DEBUG
    print_given_info(prinf, oinf);
#   endif

    struct result res = solve(prinf);
    print(res, oinf);

    free(res.s0.storage);
    free(res.s1.storage);

    return 0;
}
