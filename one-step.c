#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "asa241.h"

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
#define NUM_SAMPLESF (10.0*MAX_INFL)

double prior[NUM_SAMPLES];
double posterior[NUM_SAMPLES];
int verb = 1;


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

double variance_pdf(double *ary) {
    int i;
    double total_x = 0, total_xx = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
	double x = i/10.0;
	total_x += x*ary[i];
	total_xx += x*x*ary[i];
    }
    return total_xx - total_x*total_x;
}

double stddev_pdf(double *ary) {
    return sqrt(variance_pdf(ary));
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

void setup_beta(double *ary, double a, double b) {
    int i;
    assert(a > 0);
    assert(b > 0);
    double scale = (NUM_SAMPLES-1)/(NUM_SAMPLESF*NUM_SAMPLESF);
    for (i = 0; i < NUM_SAMPLES; i++) {
	double x = (i + 0.5)*scale;
	double p = pow(x, a-1)*pow(1 - x, b - 1);
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

void fit_minerror_norm_brute(void) {
    double best_mu = -1;
    double best_sigma = -1;
    double best_error = HUGE_VAL;
    double mu, sigma, error;
    printf("Minimum error (truncated normal) search:");
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
	if (mu - floor(mu) < 0.00001) {
	    putchar('.');
	    fflush(stdout);
	}
    }
    putchar('\n');
    printf("Minimum error fit: mu = %f, sigma = %f\n", best_mu, best_sigma);
    setup_normal(fit, best_mu, best_sigma);
}

void fit_minerror_beta_brute(void) {
    double best_a = -1;
    double best_b = -1;
    double best_error = HUGE_VAL;
    double a, b, error;
    printf("Minimum error (beta) search:");
    for (a = 0.1; a <= 100; a += 0.1) {
	for (b = 0.1; b <= 100; b += 0.1) {
	    setup_beta(fit, a, b);
	    error = calculate_error(posterior, fit);
	    if (error < best_error) {
		/* printf("Improved error to %g with mu=%f, sigma=%f\n",
		       error, mu, sigma); */
		best_a = a;
		best_b = b;
		best_error = error;
	    }
	}
	if (a - floor(a) < 0.00001) {
	    putchar('.');
	    fflush(stdout);
	}
    }
    putchar('\n');
    printf("Minimum error fit: a = %f, b = %f\n", best_a, best_b);
    setup_beta(fit, best_a, best_b);
}

void fit_moments_beta(void) {
    double mean = mean_pdf(posterior);
    double var = variance_pdf(posterior);
    double mean01 = mean/64.0;
    double var01 = var/(64.0*64.0);
    double mean1m = mean01*(1.0 - mean01);
    double diff = mean1m/var01 - 1.0;
    double a = mean01*diff;
    double b = (1.0 - mean01)*diff;
    printf("Method of moments fit: a = %f, b = %f\n", a, b);
    setup_beta(fit, a, b);
}

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

double sum_pdf(double *ary, double min, double max) {
    int start = floor(min * 10);
    int end = ceil(max * 10);
    double sum = 0;
    int i;
    assert(min <= max);
    assert(start <= end);
    for (i = start; i <= end; i++) {
	sum += ary[i];
    }
    assert(sum >= 0 && sum < 1.0001);
    return sum;
}

void fit_minsigma_brute(double posterior_mean, double posterior_stddev) {
    double dist = 2;
    double best_mu = -1;
    double best_sigma = HUGE_VAL;
    double mu, sigma;
    double min_mu = posterior_mean - dist;
    double max_mu = posterior_mean + dist;
    double min_sigma = posterior_stddev;
    double max_sigma = posterior_stddev + dist;
    
    if (min_mu < 0) {
        min_mu = 0;
    }
    if (min_sigma < 0.1) {
        min_sigma = 0.1;
    }

    for (mu = min_mu; mu <= max_mu; mu += 0.01) {
	for (sigma = min_sigma; sigma < max_sigma; sigma += 0.01) {
	    int cl_i;
	    int is_good = 1;
	    for (cl_i = 1; cl_i < 100; cl_i++) {
		double norm_lower, norm_upper;
		double cl = cl_i/100.0;
		double post_sum;
		trunc_norm_conf_interval(mu, sigma, 64, cl,
					 &norm_lower, &norm_upper);
		post_sum = sum_pdf(posterior, norm_lower, norm_upper);
		if (cl > post_sum) {
		    /* printf("sigma = %f is over-confident "
			   "at %d%% (%f vs. %f)\n",
			   sigma, cl_i, cl, post_sum); */
		    is_good = 0;
		    break;
		}
	    }
	    if (is_good)
		break;
	    else if (sigma > best_sigma)
		break;
	}
	if (sigma < best_sigma) {
	    if(verb == 1) {
            printf("For mu = %f, best safe sigma is %f\n", mu, sigma);
        }
	    best_sigma = sigma;
	    best_mu = mu;
	}
    }
    if(verb == 1) {
        printf("Minimum safe sigma is %f, for mu = %f\n", best_sigma, best_mu);
    } else {
		printf("%.2f %.2f\n", best_mu, best_sigma);
	}
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
#define PRIOR_BETA 2

int prior_type = PRIOR_UNIFORM;

#define POSTERIOR_NORMAL 1 /* truncated normal, to be precise */
#define POSTERIOR_BETA 2
#define MAX_POSTERIOR 3

#define FITNESS_MINERROR 0
#define FITNESS_MINSIGMA 1
#define FITNESS_MOMENTS 2
#define MAX_FITNESS 3

int todo[MAX_FITNESS][MAX_POSTERIOR];

#define UPDATE_EXACT 0
#define UPDATE_ATLEAST 1

int update_mode = UPDATE_EXACT;

int nSat = 0;
int k = 1;

double mu = 32;
double sigma = 18;

double a = 2;
double b = 2;

int main(int argc, char **argv) {
    int i;
    for (i = 1; i < argc; i++) {
	if (!strcmp(argv[i], "-prior") && i + 1 < argc) {
	    char *arg = argv[++i];
	    if (!strcmp(arg, "uniform")) {
		prior_type = PRIOR_UNIFORM;
	    } else if (!strcmp(arg, "normal")) {
		prior_type = PRIOR_NORMAL;
	    } else if (!strcmp(arg, "beta")) {
		prior_type = PRIOR_BETA;
	    } else {
		fprintf(stderr, "Unrecognized -prior type: "
			"should be uniform, normal, or beta\n");
		exit(1);
	    }
	} else if (!strcmp(argv[i], "-minerror") && i + 1 < argc) {
	    char *arg = argv[++i];
	    if (!strcmp(arg, "normal")) {
		todo[FITNESS_MINERROR][POSTERIOR_NORMAL] = 1;
	    } else if (!strcmp(arg, "beta")) {
		todo[FITNESS_MINERROR][POSTERIOR_BETA] = 1;
	    } else {
		fprintf(stderr, "Unrecognized -minerror type: "
			"should be normal or beta\n");
		exit(1);
	    }
	} else if (!strcmp(argv[i], "-minsigma") && i + 1 < argc) {
	    char *arg = argv[++i];
	    if (!strcmp(arg, "normal")) {
		todo[FITNESS_MINSIGMA][POSTERIOR_NORMAL] = 1;
	    } else if (!strcmp(arg, "beta")) {
		todo[FITNESS_MINSIGMA][POSTERIOR_BETA] = 1;
		assert(0); /* unimplemented */
	    } else {
		fprintf(stderr, "Unrecognized -minsigma type: "
			"should be normal or beta\n");
		exit(1);
	    }
	} else if (!strcmp(argv[i], "-moments") && i + 1 < argc) {
	    char *arg = argv[++i];
	    if (!strcmp(arg, "normal")) {
		todo[FITNESS_MOMENTS][POSTERIOR_NORMAL] = 1;
		assert(0); /* unimplemented */
	    } else if (!strcmp(arg, "beta")) {
		todo[FITNESS_MOMENTS][POSTERIOR_BETA] = 1;
	    } else {
		fprintf(stderr, "Unrecognized -moments type: "
			"should be normal or beta\n");
		exit(1);
	    }
	} else if (!strcmp(argv[i], "-mu") && i + 1 < argc) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -mu should be a number\n");
		exit(1);
	    }
	    mu = val;
	} else if (!strcmp(argv[i], "-sigma") && i + 1 < argc) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -sigma should be a number\n");
		exit(1);
	    }
	    sigma = val;
	} else if (!strcmp(argv[i], "-a") && i + 1 < argc) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -a should be a number\n");
		exit(1);
	    }
	    a = val;
	} else if (!strcmp(argv[i], "-b") && i + 1 < argc) {
	    char *arg = argv[++i];
	    char *endptr;
	    double val = strtod(arg, &endptr);
	    if (endptr == arg) {
		fprintf(stderr, "Argument to -b should be a number\n");
		exit(1);
	    }
	    b = val;
	} else if (!strcmp(argv[i], "-nSat") && i + 1 < argc) {
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
	} else if (!strcmp(argv[i], "-nSatGE") && i + 1 < argc) {
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
	} else if (!strcmp(argv[i], "-k") && i + 1 < argc) {
	    char *arg = argv[++i];
	    char *endptr;
	    long val = strtol(arg, &endptr, 0);
	    if (endptr == arg || arg < 0) {
		fprintf(stderr, "Argument to -k should be "
			"a non-negative integer\n");
		exit(1);
	    }
	    k = val;
	} else if (!strcmp(argv[i], "-verb") && i + 1 < argc) {
        char *arg = argv[++i];
        char *endptr;
        long val = strtol(arg, &endptr, 0);
        if (endptr == arg || arg < 0) {
            fprintf(stderr, "Argument to -verb should be "
                    "0 or 1\n");
            exit(1);
        }
        verb = val;
	} else {
	    fprintf(stderr, "Urecognized option `%s'\n", argv[i]);
	    exit(1);
	}
    }

    if (prior_type == PRIOR_UNIFORM) {
	setup_prior_uniform();
    } else if (prior_type == PRIOR_NORMAL) {
	setup_normal(prior, mu, sigma);
    } else if (prior_type == PRIOR_BETA) {
	setup_beta(prior, a, b);
    } else {
	assert(0);
    }
    
    if (verb == 1){
		gnuplot_data("prior.dat", prior, NUM_SAMPLES);
        printf("Prior mean is %g\n", mean_pdf(prior));
        printf("Prior stddev is %g\n", stddev_pdf(prior));
    }

    if (update_mode == UPDATE_EXACT) {
	estimate_posterior_eq_n(k, nSat);
    } else if (update_mode == UPDATE_ATLEAST) {
	estimate_posterior_ge_n(k, nSat);
    }
    
    double mean = round(mean_pdf(posterior)*100)/100;
    double stddev = round(stddev_pdf(posterior)*100)/100;
    
    if (verb == 1) {
		gnuplot_data("posterior.dat", posterior, NUM_SAMPLES);
		printf("Posterior mean is %g\n", mean_pdf(posterior));
		printf("Posterior stddev is %g\n", stddev_pdf(posterior));
	}

    if (todo[FITNESS_MOMENTS][POSTERIOR_BETA]) {
	fit_moments_beta();
	gnuplot_data("moments-beta.dat", fit, NUM_SAMPLES);
    }

    if (todo[FITNESS_MINSIGMA][POSTERIOR_NORMAL]) {
	fit_minsigma_brute(mean, stddev);
	gnuplot_data("minsigma-normal.dat", fit, NUM_SAMPLES);
    }

    if (todo[FITNESS_MINERROR][POSTERIOR_NORMAL]) {
	fit_minerror_norm_brute();
	gnuplot_data("minerror-normal.dat", fit, NUM_SAMPLES);
    }

    if (todo[FITNESS_MINERROR][POSTERIOR_BETA]) {
	fit_minerror_beta_brute();
	gnuplot_data("minerror-beta.dat", fit, NUM_SAMPLES);
    }

    return 0;
}
