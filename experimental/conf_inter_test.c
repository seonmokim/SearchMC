#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "asa241.h"

double norm_cdf(double x) {
    return 0.5*erfc(-x/sqrt(2));
}

double norm_inv_cdf(double p) {
    return r8_normal_01_cdf_inverse(p);
}

void trunc_norm_conf_interval(double mu, double sigma, double bound, double cl,
			      double *lowerp, double *upperp) {
    double denom1 = sqrt(2)*sigma;
    double integral = 0.5*(erf((bound - mu)/denom1) - erf(-mu/denom1));
    double norm = 1/integral;
    double ci_factor = norm_inv_cdf((cl/norm + 1)/2);
    double ci_sigma = ci_factor * sigma;

    double normal_upper = mu + ci_sigma;
    double normal_lower = mu - ci_sigma;

    double modif_upper =
	mu + sigma*norm_inv_cdf(cl/norm + norm_cdf(-mu/sigma));
    double modif_lower =
	mu + sigma*norm_inv_cdf(norm_cdf((bound-mu)/sigma) - cl/norm);

    double upper, lower;
    if (cl == 1.0) {
	upper = bound;
	lower = 0;
    } else if (mu > bound/2.0 && mu <= bound) {
	if (normal_upper <= bound) {
	    upper = normal_upper;
	    lower = normal_lower;
	} else {
	    upper = bound;
	    lower = modif_lower;
	}
    } else if (mu <= bound/2.0 && mu >= 0) {
	if (normal_lower >= 0) {
	    upper = normal_upper;
	    lower = normal_lower;
	} else {
	    upper = modif_upper;
	    lower = 0;
	}
    } else if (mu > bound) {
	upper = bound;
	lower = modif_lower;
    } else {
	assert(mu < 0);
	upper = modif_upper;
	lower = 0;
    }
    assert(lower >= 0 && lower <= bound);
    assert(upper >= 0 && upper <= bound);
    assert(lower <= upper);
    *lowerp = lower;
    *upperp = upper;
}

int main(int argc, char **argv) {
    double upper, lower;
    double mu, sigma;
    int i;
    if (argc == 1) {
	mu = 32;
	sigma = 18;
    } else if (argc == 3) {
	mu = atof(argv[1]);
	sigma = atof(argv[2]);
    } else {
	fprintf(stderr, "Usage: conf_inter_test <mu> <sigma>\n");
	exit(1);
    }
    for (i = 0; i <= 100; i += 5) {
	double cl = i/100.0;
	trunc_norm_conf_interval(mu, sigma, 64, cl, &lower, &upper);
	printf("%d%%: (%g, %g)\n", i, lower, upper);
    }
    return 0;
}
