#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

void setup_prior(void) {
    int i;
    /* Uniform distribution */
    for (i = 0; i < NUM_SAMPLES; i++) {
	prior[i] = 1.0;
    }
    normalize(prior, NUM_SAMPLES);
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

int main(int argc, char **argv) {
    setup_prior();
    gnuplot_data("prior.dat", prior, NUM_SAMPLES);
    printf("Prior mean is %g\n", mean_pdf(prior));
    printf("Prior stddev is %g\n", stddev_pdf(prior));
    estimate_posterior_eq_n(32, 0);
    gnuplot_data("posterior0.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat = 0) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat = 0) stddev is %g\n", stddev_pdf(posterior));
    estimate_posterior_eq_n(32, 1);
    gnuplot_data("posterior1.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat = 1) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat = 1) stddev is %g\n", stddev_pdf(posterior));
    estimate_posterior_eq_n(32, 2);
    gnuplot_data("posterior2.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat = 2) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat = 2) stddev is %g\n", stddev_pdf(posterior));
    estimate_posterior_eq_n(32, 10);
    gnuplot_data("posterior10.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat = 10) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat = 10) stddev is %g\n", stddev_pdf(posterior));

    estimate_posterior_ge_n(32, 1);
    gnuplot_data("posterior_ge1.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat >= 1) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat >= 1) stddev is %g\n", stddev_pdf(posterior));
    estimate_posterior_ge_n(32, 2);
    gnuplot_data("posterior_ge2.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat >= 2) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat >= 2) stddev is %g\n", stddev_pdf(posterior));
    estimate_posterior_ge_n(32, 10);
    gnuplot_data("posterior_ge10.dat", posterior, NUM_SAMPLES);
    printf("Posterior(nSat >= 10) mean is %g\n", mean_pdf(posterior));
    printf("Posterior(nSat >= 10) stddev is %g\n", stddev_pdf(posterior));
    return 0;
}
