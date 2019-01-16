#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

#define N_COUNT 10000001

inline double kernel(double x, double s0, double s1)
{
    double xx = x * x;
    return exp(-0.5 * (s0 * xx + s1 * xx * xx) / (1.0 + xx));
}



void calculate(struct vector_func *func, double *k, double *d)
{
    int i;
    double x = func->grid.origin;
    double norm = get_norm(func);

    for(i = 0; i < func->grid.count; i++){
        func->storage[i] *= x * x;
        x += func->grid.step;
    }

    *d = get_integral(func) / norm;

    x = func->grid.origin;
    for(i = 0; i < func->grid.count; i++){
        func->storage[i] *= x * x;
        x += func->grid.step;
    }

    *k = get_integral(func) / norm / (*d) / (*d) - 3;
}



static double get_begin(double s0, double s1)
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



int main(int argc, const char **argv)
{
    double s0;
    double s1;
    double k;
    double d;
    int i;
    double x;
    struct vector_func func;

    sscanf(argv[1], "%lf", &s0);
    sscanf(argv[2], "%lf", &s1);

    func.storage = malloc(sizeof(double) * N_COUNT);
    func.grid.origin = x = get_begin(s0, s1);
    func.grid.count = N_COUNT;
    func.grid.step = 2 * fabs(func.grid.origin) / (N_COUNT - 1);
    for(i = 0; i < N_COUNT; i++){
        func.storage[i] = kernel(x, s0, s1);
        x += func.grid.step;
    }

    calculate(&func, &k, &d);
    printf("k = %lf, d = %lf\n", k, d);

    free(func.storage);
    return 0;
}
