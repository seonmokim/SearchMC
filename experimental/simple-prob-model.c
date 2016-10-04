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

void test_choose(void) {
    int n, k;
    for (n = 0; n < 10; n++) {
	for (k = 0; k <= n; k++) {
	    printf("%g ", exp(log_choose(n, k)));
	}
	putchar('\n');
    }
}

void one_table(int bt) {
    double b = (double)bt / 10.0;
    double s = floor(pow(2.0, b) + 0.5);
    int k;
    printf("%f bits, %g solutions\n", b, s);
    for (k = 0; k <= 70; k++) {
	double rest = 1.0;
	int n;
	printf("%2d", k);
	for (n = 0; n <= 100; n++) {
	    double log_p1, log_p2, log_p, p;
	    if (n > s)
		break;
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
	    rest -= p;
	    printf(" %.4f", p);
	}
	printf(" %.4f", fabs(rest));
	putchar('\n');
    }
    putchar('\n');
}

int main(int argc, char **argv) {
    int bt;
    /* test_choose(); */

    one_table(0);
    one_table(10); one_table(16);
    one_table(20); one_table(23); one_table(26); one_table(28);
    one_table(30); one_table(32); one_table(33); one_table(34);
    for (bt = 36; bt <= 639; bt++) {
	one_table(bt);
    }
    return 0;
}
