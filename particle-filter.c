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
#define SAMPLES_PER_BIT 10
#define SAMPLES_PER_BITF ((double)SAMPLES_PER_BIT)
#define NUM_SAMPLES 500
#define NUM_SAMPLESF ((double)SAMPLES_PER_BIT*MAX_INFL)

struct map_prob {
    double influence;
    double prob;
};

struct map_prob prior[NUM_SAMPLES];
struct map_prob posterior[NUM_SAMPLES];
struct map_prob samples[NUM_SAMPLES];
double cdf[NUM_SAMPLES+1];

void normalize(struct map_prob *ary, int size) {
    int i;
    double total = 0, new_total = 0;
    for (i = 0; i < size; i++) {
        assert(ary[i].prob >= 0);
        total += ary[i].prob;
    }
    assert(total > 0);
    for (i = 0; i < size; i++) {
        ary[i].prob /= total;
        assert(ary[i].prob >= 0 && ary[i].prob <= 1.0);
        new_total += ary[i].prob;
    }
    assert(new_total >= 0.999 && new_total <= 1.001);
}

void normalize_cdf(double *ary, int size) {
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

double mean_pdf(struct map_prob *ary) {
    int i;
    double total = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
        total += ary[i].influence*ary[i].prob;
    }
    return total;
}

double variance_pdf(struct map_prob *ary) {
    int i;
    double total_x = 0, total_xx = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
        total_x += ary[i].influence*ary[i].prob;
        total_xx += ary[i].influence*ary[i].influence*ary[i].prob;
    }
    return total_xx - total_x*total_x;
}

double stddev_pdf(struct map_prob *ary) {
    return sqrt(variance_pdf(ary));
}

double mean_particle(struct map_prob *ary) {
    int i;
    double total = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
        total += ary[i].influence;
    }
    return total/NUM_SAMPLES;
}

double variance_particle(struct map_prob *ary) {
    int i;
    double total_x = 0, total_xx = 0;
    for (i = 0; i < NUM_SAMPLES; i++) {
        total_x += ary[i].influence;
        total_xx += ary[i].influence*ary[i].influence;
    }
    return total_xx/NUM_SAMPLES - total_x*total_x/NUM_SAMPLES/NUM_SAMPLES;
}

double stddev_particle(struct map_prob *ary) {
    return sqrt(variance_particle(ary));
}

void setup_prior_uniform(int nVars) {
    int i;
    for (i = 0; i < NUM_SAMPLES; i++) {
        prior[i].influence = (double)(nVars*i)/NUM_SAMPLES;
        prior[i].prob = 1.0;
    }
    normalize(prior, NUM_SAMPLES);
}

void setup_equal_weight(struct map_prob *ary) {
    int i;
    for (i = 0; i < NUM_SAMPLES; i++) {
        ary[i].prob = 1.0;
    }
    normalize(ary, NUM_SAMPLES);
}

void setup_prior_particle(char *filename) {
    FILE* file;
    if((file = fopen(filename, "r"))==NULL) {
        fprintf(stderr, "couldn't open the requested file!\n");
        exit(1);
    }
    
    int item;
    int i = 0;
    
    for (i = 0; i < NUM_SAMPLES; i++) {
        item = fscanf(file, "%lf %lf", &prior[i].influence, &prior[i].prob);
    }
    fclose(file);
    normalize(prior, NUM_SAMPLES);
}

double prob_eq_n(double bt, int k, int n) {
    double s = floor(pow(2.0, bt) + 0.5);
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
        fprintf(stderr, "bt = %g, k = %d, n = %d, s = %g\n",
                bt, k, n, s);
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

double prob_ge_n(double bt, int k, int n) {
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

void estimate_posterior_eq_n(int k, int n, struct map_prob *a, struct map_prob *b) {
    int bt;
    for (bt = 1; bt < NUM_SAMPLES; bt++) {
        double p = a[bt].prob*prob_eq_n(a[bt].influence, k, n);
        b[bt].prob = p;
        b[bt].influence = a[bt].influence;
    }
    b[0].prob = b[1].prob;
    b[0].influence = 0;
    normalize(b, NUM_SAMPLES);
}

void estimate_posterior_ge_n(int k, int n, struct map_prob *a, struct map_prob *b) {
    int bt;
    for (bt = 1; bt < NUM_SAMPLES; bt++) {
        double p = a[bt].prob*prob_ge_n(a[bt].influence, k, n);
        b[bt].prob = p;
        b[bt].influence = a[bt].influence;
    }
    b[0].prob = b[1].prob;
    b[0].influence = 0;
    normalize(b, NUM_SAMPLES);
}

int cmpfunc (const void * a, const void * b)
{
    struct map_prob *p1 = (struct map_prob *)a;
    struct map_prob *p2 = (struct map_prob *)b;
    if ( p1->influence > p2->influence) return 1;
    else if ( p1->influence < p2->influence) return -1;
    else return 0;
}

double linear_intep (double a, double b, double w1, double w2) {
    return (b*w1+a*w2)/(w1+w2);
}

void getCDF (struct map_prob *a) {
    int i;
    for (i = 0; i < NUM_SAMPLES+1; i++) {
        if(i == 0 ) {
            cdf[i] = a[i].prob/2;
        } else if (i == NUM_SAMPLES) {
            cdf[i] = a[i-1].prob/2;
        } else {
            cdf[i] = (a[i-1].prob+a[i].prob)/2;
        }
    }
    normalize_cdf(cdf, NUM_SAMPLES);
}
void sampling(struct map_prob *a, struct map_prob *b, int nVars) {
    getCDF(a);
    int i;
    for (i = 0; i < NUM_SAMPLES; i++) {
        double r1 = (double)rand() / (double)RAND_MAX;
        int index = 0;
        while (r1 > 0 && index < NUM_SAMPLES+1) {
            r1 = r1 - cdf[index];
            index++;
        }
        index--;
        double hi, lo;
        if (fabs(r1 + cdf[index]) > fabs(r1)) {
            hi = fabs(r1 + cdf[index]);
            lo = fabs(r1);
        } else {
            lo = fabs(r1 + cdf[index]);
            hi = fabs(r1);
        }
        if (index == 0) {
            b[i].influence = linear_intep(0, a[index].influence, lo, hi);
        } else if (index == NUM_SAMPLES + 1) {
             b[i].influence = linear_intep(a[index].influence, nVars, lo, hi);
        } else {
            if ( a[index-1].prob > a[index].prob) {
                b[i].influence = linear_intep(a[index-1].influence, a[index].influence, hi, lo);
            } else {
                b[i].influence = linear_intep(a[index-1].influence, a[index].influence, lo, hi);
            }
        }
    }
    
    qsort(b, NUM_SAMPLES, sizeof(struct map_prob), cmpfunc);
    //normalize(b, NUM_SAMPLES);
    setup_equal_weight(b);
}

double sum_weight(struct map_prob *ary, int min, int max) {
    double sum = 0;
    int i;
    assert(min <= max);
    for (i = min; i <= max; i++) {
        sum += ary[i].prob;
    }
    if (sum < 0 || sum >= 1.01) {
        for (i = 0; i < NUM_SAMPLES; i++) {
            printf ("%f ", ary[i].prob);
        }
        putchar('\n');
        printf ("%d %d error sum: %f\n",min, max, sum);
    }
    assert(sum >= 0 && sum < 1.01);
    return sum;
}

double sum_pdf(double *ary, int min, int max) {
    double sum = 0;
    int i;
    assert(min <= max);
    for (i = min; i <= max; i++) {
        sum += ary[i];
    }
    if (sum < 0 || sum >= 1.01) {
        for (i = 0; i < NUM_SAMPLES; i++) {
            printf ("%f ", ary[i]);
        }
        putchar('\n');
        printf ("%d %d error sum: %f\n",min, max, sum);
    }
    assert(sum >= 0 && sum < 1.01);
    return sum;
}

int getCenterIndex()
{
    int i;
    int cnt = 1;
    double max = posterior[1].prob;
    int sum = 1;
    for (i = 1 ; i < NUM_SAMPLES; i++) {
        if(max < posterior[i].prob) {
            max = posterior[i].prob;
            sum = i;
            cnt = 1;
        } else if (max == posterior[i].prob) {
            cnt++;
            sum = sum + i;
        }
    }
    //printf("sum is %g %d %d\n",max, sum ,cnt);
    int res = (int)(sum/cnt);
    return res;
}

int getIndex(double value)
{
    int res;
    int i;
    for (i = 0 ; i < NUM_SAMPLES; i++) {
        if(value < posterior[i].influence) {
            res = i;
            break;
        }
    }
    return res;
}

double getLowerBound(int cen, double confidence) {
    double res;
    double p = confidence + 0.5 * cdf[cen];
    while (p > 0 && cen > -1) {
        p = p - cdf[cen];
        cen--;
    }
    cen++;
    //printf("Lower Index is %d\n", cen);
    res = linear_intep(posterior[cen].influence, posterior[cen+1].influence, fabs(p), fabs(p + cdf[cen]));
    return res;
}

double getUpperBound(int cen, double confidence) {
    double res;
    double p = confidence + 0.5 * cdf[cen];
    while (p > 0 && cen < NUM_SAMPLES) {
        p = p - cdf[cen];
        cen++;
    }
    cen--;
    //printf("Upper Index is %d\n", cen);
    res = linear_intep(posterior[cen-1].influence, posterior[cen].influence, fabs(p), fabs(p + cdf[cen]));
    return res;
}

void getBounds(double conf) {
    
    double lb, ub;
    double lconf, uconf;
    double mean = mean_particle(posterior);
    double stddev = stddev_particle(posterior);
    int center_index = getIndex(mean);
    //int center_index = getCenterIndex();
    //double center = posterior[center_index].influence;
    //printf("Center Index is %d\n", center_index);
    //printf("sum_pdf is %g\n", sum_pdf(posterior, center_index, NUM_SAMPLES-1));
    printf("%g %g ", mean, stddev);
    //printf("%g %g ", posterior[NUM_SAMPLES/2].influence, stddev*0.67449);
    
    double sum_pdf_ub = sum_pdf(cdf, center_index, NUM_SAMPLES);
    double sum_pdf_lb = sum_pdf(cdf, 0, center_index);
    if (sum_pdf_ub < (conf / 2)) {
        uconf = sum_pdf_ub - ((1 - conf) / 2);
        lconf = conf - uconf;
        //printf("Lower conf is %g\n", lconf);
        lb = getLowerBound(center_index, lconf);
        ub = getUpperBound(center_index, uconf);
    } else if (sum_pdf_lb < (conf / 2)) {
        lconf = sum_pdf_lb - ((1 - conf) / 2);
        uconf = conf - lconf;
        //printf("Upper conf is %g\n", uconf);
        lb = getLowerBound(center_index, lconf);
        ub = getUpperBound(center_index, uconf);
    } else {
        lb = getLowerBound(center_index, conf/2);
        ub = getUpperBound(center_index, conf/2);
    }
    
    //double adjust_factor = conf / 10;
    //int interval = (int)NUM_SAMPLES*(conf+adjust_factor);
    //int lb_index = (NUM_SAMPLES-interval)/2-1;
    //int ub_index = lb_index + interval;
    //printf("%g %g\n", posterior[lb_index].influence, posterior[ub_index].influence);
    printf("%g %g\n", lb, ub);
}

void gnuplot_data(const char *fname, struct map_prob *ary, int size) {
    int i;
    int res;
    FILE *fp = fopen(fname, "w");
    assert(fp);
    for (i = 0; i < size; i++) {
        fprintf(fp, "%g %g\n", ary[i].influence, ary[i].prob);
    }
    res = fclose(fp);
    assert(res == 0);
}

#define PRIOR_UNIFORM 0
#define PRIOR_PARTICLE 1

int prior_type = PRIOR_UNIFORM;

#define UPDATE_EXACT 0
#define UPDATE_ATLEAST 1

int update_mode = UPDATE_EXACT;

int nSat = 0;
int k = 1;

double center = 32;
double cl = 0;
int nVars = 64;

int verb = 1;

char *prior_file;
char *pid;

int main(int argc, char **argv) {
    int i;
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-prior") && i + 1 < argc) {
            char *arg = argv[++i];
            if (!strcmp(arg, "uniform")) {
                prior_type = PRIOR_UNIFORM;
            } else if (!strcmp(arg, "particle")) {
                prior_type = PRIOR_PARTICLE;
            } else {
                fprintf(stderr, "Unrecognized -prior type: "
                        "should be uniform, normal, beta or particle\n");
                exit(1);
            }
        } else if (!strcmp(argv[i], "-priorparticle") && i + 1 < argc) {
            prior_file = argv[++i];
        } else if (!strcmp(argv[i], "-cl") && i + 1 < argc) {
            char *arg = argv[++i];
            char *endptr;
            double val = strtod(arg, &endptr);
            if (endptr == arg) {
                fprintf(stderr, "Argument to -cl should be a number\n");
                exit(1);
            }
            cl = val;
        } else if (!strcmp(argv[i], "-center") && i + 1 < argc) {
            char *arg = argv[++i];
            char *endptr;
            double val = strtod(arg, &endptr);
            if (endptr == arg) {
                fprintf(stderr, "Argument to -center should be a number\n");
                exit(1);
            }
            center = val;
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
        } else if (!strcmp(argv[i], "-nVars") && i + 1 < argc) {
            char *arg = argv[++i];
            char *endptr;
            long val = strtol(arg, &endptr, 0);
            if (endptr == arg || arg < 0) {
                fprintf(stderr, "Argument to -nVars should be "
                        "a non-negative integer\n");
                exit(1);
            }
            nVars = val;
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
        } else if (!strcmp(argv[i], "-pid") && i + 1 < argc) {
            pid = argv[++i];
        } else {
            fprintf(stderr, "Urecognized option `%s'\n", argv[i]);
            exit(1);
        }
    }
    
    if (prior_type == PRIOR_UNIFORM) {
        setup_prior_uniform(nVars);
    } else if (prior_type == PRIOR_PARTICLE) {
        setup_prior_particle(prior_file);
    } else {
        assert(0);
    }
    
    //gnuplot_data("prior.dat", prior, NUM_SAMPLES);
    
    if (verb == 1){
        printf("Prior mean is %g\n", mean_pdf(prior));
        printf("Prior stddev is %g\n", stddev_pdf(prior));
    }
    
    if (update_mode == UPDATE_EXACT) {
        estimate_posterior_eq_n(k, nSat, prior, prior);
    } else if (update_mode == UPDATE_ATLEAST) {
        estimate_posterior_ge_n(k, nSat, prior, prior);
    }
    sampling(prior, posterior, nVars);
    if (verb == 1){
        printf("Posterior mean is %g\n", mean_pdf(posterior));
        printf("Posterior stddev is %g\n", stddev_pdf(posterior));
    }
    getCDF(posterior);
    char posterior_filename[80];
    strcpy(posterior_filename, "posterior_");
    strcat(posterior_filename, pid);
    strcat(posterior_filename, ".dat");
    
    
    gnuplot_data(posterior_filename, posterior, NUM_SAMPLES);

    getBounds(cl);

    return 0;
}
