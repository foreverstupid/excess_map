#include "vector.h"

/*
 * Weight in a quadratic integrate formula
 */
inline static double weight(int i, int n, double step)
{
    return i == 0 || i == n - 1 ? step / 2 : step;
}



double get_norm(const struct vector_func *f)
{
    int i;
    double res = 0.0;

    for(i = 0; i < f->grid.count; i++){
        res += fabs(f->storage[i]) * weight(i, f->grid.count, f->grid.step);
    }

    return res;
}



double get_integral(const struct vector_func *f)
{
    int i;
    double res = 0.0;

    for(i = 0; i < f->grid.count; i++){
        res += f->storage[i] * weight(i, f->grid.count, f->grid.step);
    }

    return res;
}
