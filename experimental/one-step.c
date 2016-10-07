#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double log_fact(double n) {
    return lgamma(n + 1);
}

/* This algorithm is constant time, but it doesn't work well when n >>
   k, because then log_fact(n) and log_fact(n - k) are too close
   together */
double log_choose_gamma(double n, double k) {
    return (log_fact(n) - log_fact(n - k)) - log_fact(k);
}

/* When n >> k, n!/(n-k)! is almost n^k, which suggests the following
   refinement of log_choose_gamma. However some more thought is
   probably needed about how to choose between the two formulas. */
double log_choose_gamma2(double n, double k) {
    if (n / k > 1000000000) {
	return log(n)*k - log_fact(k);
    } else {
	return (log_fact(n) - log_fact(n - k)) - log_fact(k);
    }
}

/* This algorithm takes O(k) steps, but as long as k is small, it
   gives good results even for very large n.
*/
double log_choose_small_k(double n, int k) {
    double l = 0;
    int i;
    for (i = 0; i < k; i++) {
	l += log(n - i);
	l -= log(k - i);
    }
    return l;
}

#define log_choose log_choose_small_k

#define MAX_INFL 64
#define NUM_SAMPLES (10*MAX_INFL)

double prior[NUM_SAMPLES];
double posterior[NUM_SAMPLES];

void normalize(double *ary, int size) {
    int i;
    double total = 0, new_total = 0;
    for (i = 0; i < size; i++) {
	assert(ary[i] >= 0);
	total += ary[i];
    }
    assert(total > 0);
    for (i = 0; i < size; i++) {
	ary[i] /= total;
	assert(ary[i] >= 0 && ary[i] <= 1.0);
	new_total += ary[i];
    }
    assert(new_total >= 0.999 && new_total <= 1.001);
}

double mean_pdf(double *ary) {
    int i;
    double total = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
	double x = i/10.0;
	total += x*ary[i];
    }
    return total;
}

double stddev_pdf(double *ary) {
    int i;
    double total_x = 0, total_xx = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
	double x = i/10.0;
	total_x += x*ary[i];
	total_xx += x*x*ary[i];
    }
    return sqrt(total_xx - total_x*total_x);
}

double prob_eq_n(int bt, int k, int n) {
    double b = (double)bt / 10.0;
    double s = floor(pow(2.0, b) + 0.5);
    double log_p1, log_p2, log_p, p;
    if (n > s)
	return 0;
    log_p1 = -k * log(2);
    if (k == 0) {
	p = (s == n) ? 1.0 : 0.0;
    } else {
	log_p2 = log1p(-pow(2.0, -k));
	log_p = log_choose(s, n) + log_p1 * n + log_p2 * (s - n);
	p = exp(log_p);
    }
    if (p < 0.0) {
	fprintf(stderr, "Bug: negative probability\n");
	exit(1);
    } else if (p > 1.0) {
	fprintf(stderr, "Bug: probability more than 1.0\n");
	fprintf(stderr, "b = %g, k = %d, n = %d, s = %g\n",
		b, k, n, s);
	fprintf(stderr, "p = %g\n", p);
	if (n == 1) {
	    fprintf(stderr, "log_choose = %g, s/b %g\n",
		    log_choose(s, n), log(s));
	}
	exit(1);
    }
    assert(p >= 0.0 && p <= 1.0);
    return p;
}

double prob_ge_n(int bt, int k, int n) {
    int n2;
    double prob = 0.0;
    for (n2 = 0; n2 < n; n2++) {
	prob += prob_eq_n(bt, k, n2);
    }
    assert(prob >= -0.0001 && prob <= 1.0001);
    if (prob < 0)
	prob = 0.0;
    if (prob > 1)
	prob = 1.0;
    return 1.0 - prob;
}

void setup_prior_uniform(void) {
    int i;
    for (i = 0; i < NUM_SAMPLES; i++) {
	prior[i] = 1.0;
    }
    normalize(prior, NUM_SAMPLES);
}

void setup_normal(double *ary, double mu, double sigma) {
    int i;
    double denom = 2 * sigma * sigma;
    assert(denom > 0.0);
    for (i = 0; i < NUM_SAMPLES; i++) {
	double x = i/10.0;
	double diff = x - mu;
	double p = exp(-(diff*diff)/denom);
	assert(p >= 0.0 && p <= 1.0);
	ary[i] = p;
    }
    normalize(ary, NUM_SAMPLES);
}

void estimate_posterior_eq_n(int k, int n) {
    int bt;
    for (bt = 0; bt < NUM_SAMPLES; bt++) {
	double p = prior[bt]*prob_eq_n(bt, k, n);
	posterior[bt] = p;
    }
    normalize(posterior, NUM_SAMPLES);
}

void estimate_posterior_ge_n(int k, int n) {
    int bt;
    for (bt = 0; bt < NUM_SAMPLES; bt++) {
	double p = prior[bt]*prob_ge_n(bt, k, n);
	posterior[bt] = p;
    }
    normalize(posterior, NUM_SAMPLES);
}

double fit[NUM_SAMPLES];

double calculate_error(double *a1, double *a2) {
    int i;
    double error = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
	double diff = a1[i] - a2[i];
	error += diff*diff;
    }
    return error;
}

void fit_minerror_brute(void) {
    double best_mu = -1;
    double best_sigma = -1;
    double best_error = HUGE_VAL;
    double mu, sigma, error;
    for (mu = 0.0; mu <= 64; mu += 0.1) {
	/* mu values outside [0, 64] would also be possible to consider,
	   but the combination of a out-of-bounds mu and a small sigma
	   can cause all the in-bounds probabilities to underflow to 0,
	   which needs to be avoided. */
	for (sigma = 0.1; sigma < 100; sigma += 0.1) {
	    setup_normal(fit, mu, sigma);
	    error = calculate_error(posterior, fit);
	    if (error < best_error) {
		/* printf("Improved error to %g with mu=%f, sigma=%f\n",
		       error, mu, sigma); */
		best_mu = mu;
		best_sigma = sigma;
		best_error = error;
	    }
	}
    }
    printf("Best fit: mu = %f, sigma = %f\n", best_mu, best_sigma);
    setup_normal(fit, best_mu, best_sigma);
}

void gnuplot_data(const char *fname, double *ary, int size) {
    int i;
    int res;
    FILE *fp = fopen(fname, "w");
    assert(fp);
    for (i = 0; i < size; i++) {
	fprintf(fp, "%g %g\n", i/10.0, ary[i]);
    }
    res = fclose(fp);
    assert(res == 0);
}

#define PRIOR_UNIFORM 0
#define PRIOR_NORMAL 1 /* truncated normal, to be precise */

int prior_type = PRIOR_UNIFORM;

#define UPDATE_EXACT 0
#define UPDATE_ATLEAST 1

int update_mode = UPDATE_EXACT;

int nSat = 0;
int k = 1;

double mu = 32;
double sigma = 18;

int main(int argc, char **argv) {
    int i;
    for (i = 1; i < argc; i++) {
	if (!strcmp(argv[i], "-prior") && i + 1 < argc) {
	    char *arg = argv[++i];
	    if (!strcmp(arg, "uniform")) {
		prior_type = PRIOR_UNIFORM;
	    } else if (!strcmp(arg, "normal")) {
		prior_type = PRIOR_NORMAL;
	    } else {
		fprintf(stderr, "Unrecognized -prior type: "
			"should be uniform or normal\n");
		exit(1);
	    }
	} else if (!strcmp(argv[i], "-mu")) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -mu should be a number\n");
		exit(1);
	    }
	    mu = val;
	} else if (!strcmp(argv[i], "-sigma")) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -sigma should be a number\n");
		exit(1);
	    }
	    sigma = val;
	} else if (!strcmp(argv[i], "-nSat")) {
	    char *arg = argv[++i];
	    char *endptr;
	    long val = strtol(arg, &endptr, 0);
	    if (endptr == arg || arg < 0) {
		fprintf(stderr, "Argument to -nSat should be "
			"a non-negative integer\n");
		exit(1);
	    }
	    nSat = val;
	    update_mode = UPDATE_EXACT;
	} else if (!strcmp(argv[i], "-nSatGE")) {
	    char *arg = argv[++i];
	    char *endptr;
	    long val = strtol(arg, &endptr, 0);
	    if (endptr == arg || arg < 0) {
		fprintf(stderr, "Argument to -nSatGE should be "
			"a non-negative integer\n");
		exit(1);
	    }
	    nSat = val;
	    update_mode = UPDATE_ATLEAST;
	} else if (!strcmp(argv[i], "-k")) {
	    char *arg = argv[++i];
	    char *endptr;
	    long val = strtol(arg, &endptr, 0);
	    if (endptr == arg || arg < 0) {
		fprintf(stderr, "Argument to -k should be "
			"a non-negative integer\n");
		exit(1);
	    }
	    k = val;
	} else {
	    fprintf(stderr, "Urecognized option `%s'\n", argv[i]);
	    exit(1);
	}
    }

    if (prior_type == PRIOR_UNIFORM) {
	setup_prior_uniform();
    } else if (prior_type == PRIOR_NORMAL) {
	setup_normal(prior, mu, sigma);
    } else {
	assert(0);
    }

    gnuplot_data("prior.dat", prior, NUM_SAMPLES);
    printf("Prior mean is %g\n", mean_pdf(prior));
    printf("Prior stddev is %g\n", stddev_pdf(prior));

    if (update_mode == UPDATE_EXACT) {
	estimate_posterior_eq_n(k, nSat);
    } else if (update_mode == UPDATE_ATLEAST) {
	estimate_posterior_ge_n(k, nSat);
    }

    gnuplot_data("posterior.dat", posterior, NUM_SAMPLES);
    printf("Posterior mean is %g\n", mean_pdf(posterior));
    printf("Posterior stddev is %g\n", stddev_pdf(posterior));

    fit_minerror_brute();
    gnuplot_data("fit.dat", fit, NUM_SAMPLES);

    return 0;
}
