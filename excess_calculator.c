#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

#define N_COUNT 10000001

typedef double (*Kernel)(double, double, double);

double exp_rat_kernel(double x, double s0, double s1)
{
    double xx = x * x;
    return exp(-0.5 * (s0 * xx + s1 * xx * xx) / (1.0 + xx));
}

double roughgarden_kernel(double x, double s, double g)
{
    return exp(-pow(fabs(x / s), g));    
}

struct kern_type {
    char type;
    Kernel kernel;
} kern_types[] = {
    { 'k', &exp_rat_kernel },
    { 'r', &roughgarden_kernel }
};



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



static double get_begin(Kernel kernel, double a, double b)
{
    double x = 0.0;
    double step = 1e-5;

    while(fabs(kernel(x, a, b)) > 1e-12){
        x -= step;
    }

#   ifdef DEBUG
    printf("Origin: %lf\n", x);
#   endif

    return x;
}



Kernel get_kernel(char kern_type)
{
    unsigned int i;
    unsigned int kern_count = sizeof(kern_types) / sizeof(struct kern_type);

    for(i = 0; i < kern_count; i++){
        if(kern_types[i].type == kern_type){
            return kern_types[i].kernel;
        }
    }

    return NULL;
}



int main(int argc, const char **argv)
{
    double a;
    double b;
    char kern;
    double k;
    double d;
    int i;
    double x;
    struct vector_func func;
    Kernel kernel; 

    kern = argv[1][0];
    kernel = get_kernel(kern);
    if(kernel == NULL){
        fprintf(stderr, "### Unknown kernel type!\n");
        return 1;
    }

    sscanf(argv[2], "%lf", &a);
    sscanf(argv[3], "%lf", &b);

    func.storage = malloc(sizeof(double) * N_COUNT);
    func.grid.origin = x = get_begin(kernel, a, b);
    func.grid.count = N_COUNT;
    func.grid.step = 2 * fabs(func.grid.origin) / (N_COUNT - 1);
    for(i = 0; i < N_COUNT; i++){
        func.storage[i] = kernel(x, a, b);
        x += func.grid.step;
    }

    calculate(&func, &k, &d);
    printf("k = %lf, d = %lf\n", k, d);

    free(func.storage);
    return 0;
}
