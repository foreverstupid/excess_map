#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include "kurtic.h"
#include "vector.h"
#include "solver.h"

typedef int (*func)(void);                 /* type of test function */

/*
 * Holds info about test function
 */
struct func_info{
    func function;
    const char *name;
};


/*
 * Holds info about success of test run
 */
enum { passed, failed };



/*======================================================================*/
/*                                ASSERTIONS                            */
/*======================================================================*/
/*
 * Check two double values equality with given precision. In case of fail
 * writes info to the standart output
 */
int assert_double(double x, double y, double eps, const char *msg)
{
    int flag = fabs(x - y) < eps;

    if(!flag){
        printf("\n%s\n    Expected: %lf\n    Actual:   %lf\n", msg, x, y);
    }

    return flag;
}

/*
 * Checks boolean value. In case of fail writes info to the standart
 * output
 */
int assert_bool(int flag, const char *msg)
{
    if(!flag){
        printf("\n%s\n    Expected: true\n    Actual:   false\n", msg);
    }

    return flag;
}

/*
 * Checks two integer values equality. In case of fail writes info to
 * the standart output
 */
int assert_int(int x, int y, const char *msg)
{
    int flag = x == y;

    if(!flag){
        printf("\n%s\n    Expected: %d\n    Actual:   %d\n", msg, x, y);
    }

    return flag;
}





/*======================================================================*/
/*                              TEST FUNCTIONS                          */
/*======================================================================*/
/*
 * Tests integral calculation
 */
int test_integral()
{
    int count = 1000;
    double origin = 0.0;
    double last = 1.0;

    double expected = M_E - 1.0;
    double actual;
    double eps = 1e-6;

    int i;
    double x;
    int result;
    struct vector_func f;

    f.storage = malloc(sizeof(double) * count);
    f.grid.count = count;
    f.grid.origin = origin;
    f.grid.step = (last - origin) / (count - 1);

    x = f.grid.origin;
    for(i = 0; i < f.grid.count; i++){
        f.storage[i] = exp(x);
        x += f.grid.step;
    }

    actual = get_integral(&f);
    result =
        assert_double(expected, actual, eps, "Integral value")
        ? passed
        : failed;

    free(f.storage);
    return result;
}



/*
 * Tests L1 norm calculation
 */
int test_norm()
{
    int count = 10000;
    double origin = 0.0;
    double last = 2 * M_PI;

    double expected = 4;
    double actual;
    double eps = 1e-6;

    int i;
    double x;
    int result;
    struct vector_func f;

    f.storage = malloc(sizeof(double) * count);
    f.grid.count = count;
    f.grid.origin = origin;
    f.grid.step = (last - origin) / (count - 1);

    x = f.grid.origin;
    for(i = 0; i < f.grid.count; i++){
        f.storage[i] = sin(x);
        x += f.grid.step;
    }

    actual = get_norm(&f);
    result =
        assert_double(expected, actual, eps, "Norm value")
        ? passed
        : failed;

    free(f.storage);
    return result;
}



/*
 * Tests solver
 */
int test_solver()
{
    int flag;
    struct problem_info pinf;
    struct result res;

    pinf.k_grid.origin = 0.0;
    pinf.k_grid.step = 0.0;
    pinf.k_grid.count = 1;

    pinf.d_grid.origin = M_PI;
    pinf.d_grid.step = 0.0;
    pinf.d_grid.count = 1;

    pinf.space_grid.origin = -10.0;
    pinf.space_grid.count = 10000;
    pinf.space_grid.step =
        (10.0 - pinf.space_grid.origin) / (pinf.space_grid.count - 1);

    res = solve(pinf);

    flag =
        assert_int(res.s0.grid.count, 1, "S0 count") &&
        assert_int(res.s1.grid.count, 1, "S1 count") &&
        assert_double(1.0 / M_PI, res.s0.storage[0], 1e-3, "s0") &&
        assert_double(1.0 / M_PI, res.s1.storage[0], 1e-3, "s1");

    free(res.s0.storage);
    free(res.s1.storage);

    return flag ? passed : failed;
}





/*======================================================================*/
/*                                   MAIN                               */
/*======================================================================*/
/*
 * Runs all tests
 */
int main(int argc, const char **argv)
{
    struct func_info test_funcs[] = {
        { &test_integral, "test_integral" },
        { &test_norm, "test_norm" },
        { &test_solver, "test_solver" }
    };

    unsigned int i;
    for(i = 0; i < sizeof(test_funcs) / sizeof(struct func_info); i++){
        printf("%s: ", test_funcs[i].name);
        if(test_funcs[i].function() == passed){
            printf("passed\n");
        }else{
            printf("FAILED!\n\n");
        }
    }

    return 0;
}
