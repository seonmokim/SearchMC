#include "asa241.h"
#include <stdio.h>

double inv_cdf(double p) {
    return r8_normal_01_cdf_inverse(p);
}

int main(int argc, char **argv) {
    printf("inv_cdf(%g) = %g\n", 0.5, inv_cdf(0.5));
    printf("inv_cdf(%g) = %g\n", 0.975, inv_cdf(0.975));
    return 0;
}
