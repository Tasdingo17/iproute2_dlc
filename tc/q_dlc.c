// SPDX-License-Identifier: GPL-2.0-only
/*
 * Author:    Stephen Hemminger <shemminger@linux-foundation.org>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <string.h>
#include <errno.h>

#include "utils.h"
#include "tc_util.h"
#include "tc_common.h"
#include "q_dlc_spec.h"

static void explain(void)
{
    fprintf(stderr,
        "Usage: ... dlc_qdisc [ limit PACKETS ]\n"
        "                 [ delay TIME JITTER [JITTER_STEPS] ]\n"
        "                 [ distribution {uniform|normal|pareto|paretonormal} ]\n"
        "                 [ loss PERCENT ]\n"
        "                 [ mu PERCENT ]\n"
        "                 [ mean_burst_len NUM ]\n"
        "                 [ mean_good_burst_len NUM ]\n"
        "                 [ rate RATE ]\n"
    );
}

static void explain1(const char *arg)
{
    fprintf(stderr, "Illegal \"%s\"\n", arg);
}

/* Upper bound on size of distribution
 *  really (TCA_BUF_MAX - other headers) / sizeof (__s16)
 */
#define MAX_DIST    (16*1024)


/* Print values only if they are non-zero */
/*
static void __attribute__((format(printf, 2, 0)))
__print_int_opt(const char *label_json, const char *label_fp, int val)
{
    print_int(PRINT_JSON, label_json, NULL, val);
    if (val != 0)
        print_int(PRINT_FP, NULL, label_fp, val);
}
#define PRINT_INT_OPT(label, val)            \
    __print_int_opt(label, " " label " %d", (val))
*/

/* Time print prints normally with varying units, but for JSON prints
 * in seconds (1ms vs 0.001).
 */
/*
static void __attribute__((format(printf, 2, 0)))
__print_time64(const char *label_json, const char *label_fp, __u64 val)
{
    SPRINT_BUF(b1);

    print_string(PRINT_FP, NULL, label_fp, sprint_time64(val, b1));
    print_float(PRINT_JSON, label_json, NULL, val / 1000000000.);
}
#define __PRINT_TIME64(label_json, label_fp, val)    \
    __print_time64(label_json, label_fp " %s", (val))
#define PRINT_TIME64(label, val) __PRINT_TIME64(label, " " label, (val))
*/

/* Percent print prints normally in percentage points, but for JSON prints
 * an absolute value (1% vs 0.01).
 */
static void __attribute__((format(printf, 2, 0)))
__print_percent(const char *label_json, const char *label_fp, __u32 per)
{
    print_float(PRINT_FP, NULL, label_fp, (100. * per) / UINT32_MAX);
    print_float(PRINT_JSON, label_json, NULL, (1. * per) / UINT32_MAX);
}
#define __PRINT_PERCENT(label_json, label_fp, per)        \
    __print_percent(label_json, label_fp " %g%%", (per))
#define PRINT_PERCENT(label, per) __PRINT_PERCENT(label, " " label, (per))

/* scaled value used to percent of maximum. */
static void set_percent(__u32 *percent, double per)
{
    *percent = rint(per * UINT32_MAX);
}

static int get_percent(__u32 *percent, const char *str)
{
    double per;

    if (parse_percent(&per, str))
        return -1;

    set_percent(percent, per);
    return 0;
}

static void print_transition_probs(__u16 probs[3])
{
    if (!is_json_context())
        return;
    char res[64];
    sprintf(res, "[%d, %d, %d]", probs[0], probs[1], probs[2]);
    print_string(PRINT_JSON, "trans_probs", NULL, res);
}

/*
 * Simplistic file parser for distribution data.
 * Format is:
 *    # comment line(s)
 *    data0 data1 ...
 */
static int get_distribution(const char *type, __s16 *data, int maxdata)
{
    FILE *f;
    int n;
    long x;
    size_t len;
    char *line = NULL;
    char name[128];

    snprintf(name, sizeof(name), "%s/%s.dist", get_tc_lib(), type);
    f = fopen(name, "r");
    if (f == NULL) {
        fprintf(stderr, "No distribution data for %s (%s: %s)\n",
            type, name, strerror(errno));
        return -1;
    }

    n = 0;
    while (getline(&line, &len, f) != -1) {
        char *p, *endp;

        if (*line == '\n' || *line == '#')
            continue;

        for (p = line; ; p = endp) {
            x = strtol(p, &endp, 0);
            if (endp == p)
                break;

            if (n >= maxdata) {
                fprintf(stderr, "%s: too much data\n",
                    name);
                n = -1;
                goto error;
            }
            data[n++] = x;
        }
    }
 error:
    free(line);
    fclose(f);
    return n;
}

#define NEXT_IS_NUMBER() (NEXT_ARG_OK() && isdigit(argv[1][0]))
#define NEXT_IS_SIGNED_NUMBER() \
    (NEXT_ARG_OK() && (isdigit(argv[1][0]) || argv[1][0] == '-'))


static __u32 _adjust_jitter_steps(__s64 latency, __s64 jitter, __u32 old_jitter_steps, double loss_perc, double mu_perc){
    // pi1 + pi2 + pi3 = 1; markov chain states probabilities
    double pi2 = mu_perc * (1 - loss_perc);
    double pi1 = (1 - loss_perc) * (1 - mu_perc);

    __u32 new_steps = old_jitter_steps;
    __s64 step = jitter / new_steps;
    double q_mean_delay = (latency - pi1 * (latency - step)) / pi2;
    while ((__s64) q_mean_delay > (latency +jitter)){
        new_steps += 1;
        step = jitter / new_steps;
        q_mean_delay = (latency - pi1 * (latency - step)) / pi2;
    }
    if (new_steps != old_jitter_steps){
        fprintf(stderr, "Warning: increase jitter_step=%d to adjust mean_delay in queue_state", new_steps);
    }

    return new_steps;
}

struct MM1CalcParams {
    __u32 j_steps;
    double offset;
};


static double newton_method(double (*f)(double, const struct MM1CalcParams*),
                     double (*df)(double, const struct MM1CalcParams*),
                     const struct MM1CalcParams* params, 
                     double initial_guess, 
                     double tolerance, 
                     int max_iter
){
    double x = initial_guess;
    double fx = 0;
    for (int i = 0; i < max_iter; i++) {
        fx = f(x, params);
        if (fabs(fx) < tolerance)
            break;
        double dfx = df(x, params);
        x = x - fx/dfx;
    }
    if (!(fabs(fx) < tolerance)){
        fprintf(stderr, "Warning: newton method for mm1k_rho didn't converge");
    }
    return x;
};

static double func(double x, const struct MM1CalcParams* params) {
    double x_pw = pow(x, params->j_steps + 1);
    return x / (1-x) - (params->j_steps + 1) * x_pw / (1 - x_pw) - params->offset;
};

static double derivative(double x, const struct MM1CalcParams* params) {
    double x_pw = pow(x, params->j_steps + 1);
    double term1 = 1/(1-x) + x/pow(1-x, 2);
    double term2 = pow(params->j_steps+1, 2) * pow(x, params->j_steps)/(1 - x_pw);
    double term3 = pow(params->j_steps+1, 2) * pow(x, 2*params->j_steps + 1)/pow(1 - x_pw, 2);
    return term1 - term2 - term3;
};

static __u32 _calc_mm1k_rho(__s64 latency, __s64 jitter, __u32 jitter_steps, double loss_perc, double mu_perc){
    // pi1 + pi2 + pi3 = 1; markov chain states probabilities
    double pi2 = mu_perc * (1 - loss_perc);
    double pi1 = (1 - loss_perc) * (1 - mu_perc);
    __s64 step = jitter / jitter_steps;
    double q_mean_delay = (latency - pi1 * (latency - step)) / pi2;
    struct MM1CalcParams params = {
        .j_steps = jitter_steps,
        .offset = (q_mean_delay - latency) / step
    };

    double init_guess = 0.5;
    double tolerance = 0.000001;
    int max_iter = 10000;
    return (__u32) newton_method(func, derivative, &params, init_guess, tolerance, max_iter);
}

static int dlc_parse_opt(const struct qdisc_util *qu, int argc, char **argv,
               struct nlmsghdr *n, const char *dev)
{
    int dist_size = 0;
    struct rtattr *tail;
    struct tc_dlc_qopt opt = { .limit = 1000, .jitter_steps = 16 };
    __s16 *dist_data = NULL;
    int present[__TCA_DLC_MAX] = {};
    __s64 latency64 = 0;
    __s64 jitter64 = 0;
    __u64 rate64 = 0;

    double loss_perc = 0, mu_perc = 0;  // for internal calculations, vals in [0, 1]

    for ( ; argc > 0; --argc, ++argv) {
        if (matches(*argv, "limit") == 0) {
            NEXT_ARG();
            if (get_size(&opt.limit, *argv)) {
                explain1("limit");
                return -1;
            }
        } else if (matches(*argv, "latency") == 0 || matches(*argv, "delay") == 0) {
            NEXT_ARG();

            /* Old latency value in opt is no longer used. */
            present[TCA_DLC_LATENCY64] = 1;

            if (get_time64(&latency64, *argv)) {
                explain1("latency");
                return -1;
            }

            if (!NEXT_IS_NUMBER()) {
                explain1("jitter");
                return -1;
            } 
            NEXT_ARG();
            present[TCA_DLC_JITTER64] = 1;
            if (get_time64(&jitter64, *argv)) {
                explain1("latency");
                return -1;
            }

            if (NEXT_IS_NUMBER()) {
                NEXT_ARG();
                if (get_u32(&opt.jitter_steps, *argv, 0)) {
                    explain1("jitter_steps");
                    return -1;
                }
            }
            
        } else if (matches(*argv, "loss") == 0 || matches(*argv, "drop") == 0) {
            NEXT_ARG();
            if (get_percent(&opt.loss, *argv)) {
                explain1("loss percent");
                return -1;
            }
            if (parse_percent(&loss_perc, *argv)){
                return -1;
            }
            
        } else if (matches(*argv, "mu") == 0) {
            NEXT_ARG();
            if (get_percent(&opt.mu, *argv)) {
                explain1("mu");
                return -1;
            }
            if (parse_percent(&mu_perc, *argv)){
                return -1;
            }
            
        } else if (matches(*argv, "mean_burst_len") == 0) {
            NEXT_ARG();
            if (get_u32(&opt.mean_burst_len, *argv, 0)) {
                explain1("mean_burst_len");
                return -1;
            }
        } else if (matches(*argv, "mean_good_burst_len") == 0) {
            NEXT_ARG();
            if (get_u32(&opt.mean_good_burst_len, *argv, 0)) {
                explain1("mean_good_burst_len");
                return -1;
            }
        } else if (matches(*argv, "distribution") == 0) {
            NEXT_ARG();
            dist_data = calloc(MAX_DIST, sizeof(dist_data[0]));
            if (dist_data == NULL)
                return -1;

            dist_size = get_distribution(*argv, dist_data, MAX_DIST);
            if (dist_size <= 0) {
                free(dist_data);
                return -1;
            }
        } else if (matches(*argv, "rate") == 0) {
            present[TCA_DLC_RATE64] = 1;
            NEXT_ARG();
            if (strchr(*argv, '%')) {
                if (get_percent_rate64(&rate64, *argv, dev)) {
                    explain1("rate");
                    return -1;
                }
            } else if (get_rate64(&rate64, *argv)) {
                explain1("rate");
                return -1;
            }
        } else {
            if (strcmp(*argv, "help") != 0)
                fprintf(stderr, "What is \"%s\"?\n", *argv);
            explain();
            return -1;
        }
    }

    tail = NLMSG_TAIL(n);

    if (dist_data && (latency64 == 0 || jitter64 == 0)) {
        fprintf(stderr, "distribution specified but no latency and jitter values\n");
        explain();
        return -1;
    }

    opt.jitter_steps = _adjust_jitter_steps(latency64, jitter64, opt.jitter_steps, loss_perc, mu_perc);
    opt.mm1_rho = _calc_mm1k_rho(latency64, jitter64, opt.jitter_steps, loss_perc, mu_perc);

    if (addattr_l(n, 1024, TCA_OPTIONS, &opt, sizeof(opt)) < 0)
        return -1;

    if (present[TCA_DLC_LATENCY64] &&
        addattr_l(n, 1024, TCA_DLC_LATENCY64, &latency64, sizeof(latency64)) < 0)
        return -1;

    if (present[TCA_DLC_JITTER64] &&
        addattr_l(n, 1024, TCA_DLC_JITTER64, &jitter64, sizeof(jitter64)) < 0)
        return -1;
    
    if (rate64 >= (1ULL << 32)) {
        if (addattr_l(n, 1024, TCA_DLC_RATE64, &rate64, sizeof(rate64)) < 0){
            return -1;
        }
        opt.rate = ~0U;
    } else {
        opt.rate = rate64;
    }

    if (dist_data) {
        if (addattr_l(n, MAX_DIST * sizeof(dist_data[0]),
                  TCA_DLC_DELAY_DIST,
                  dist_data, dist_size * sizeof(dist_data[0])) < 0)
            return -1;
        free(dist_data);
    }

    tail->rta_len = (void *) NLMSG_TAIL(n) - (void *) tail;
    return 0;
}

static int dlc_print_opt(const struct qdisc_util *qu, FILE *f, struct rtattr *opt)
{
    struct tc_dlc_qopt qopt;
    int len;
    __u64 rate64 = 0;
    __u64 latency64 = 0;
    __u64 jitter64 = 0;
    struct tc_dlc_model* model = NULL;
    struct tc_dlc_simple_state* simple_state = NULL;
    struct tc_dlc_queue_state* queue_state = NULL;
    struct tc_dlc_loss_state* loss_state = NULL;

    if (opt == NULL)
        return 0;

    len = RTA_PAYLOAD(opt) - sizeof(qopt);
    if (len < 0) {
        fprintf(stderr, "options size error\n");
        return -1;
    }
    memcpy(&qopt, RTA_DATA(opt), sizeof(qopt));

    if (len > 0) {
        struct rtattr *tb[TCA_DLC_MAX+1];

        parse_rtattr(tb, TCA_DLC_MAX, RTA_DATA(opt) + sizeof(qopt), len);

        if (tb[TCA_DLC_RATE64]) {
            if (RTA_PAYLOAD(tb[TCA_DLC_RATE64]) < sizeof(rate64))
                return -1;
            rate64 = rta_getattr_u64(tb[TCA_DLC_RATE64]);
        }
        if (tb[TCA_DLC_LATENCY64]) {
            if (RTA_PAYLOAD(tb[TCA_DLC_LATENCY64]) < sizeof(latency64))
                return -1;
            latency64 = rta_getattr_u64(tb[TCA_DLC_LATENCY64]);

        }
        if (tb[TCA_DLC_JITTER64]) {
            if (RTA_PAYLOAD(tb[TCA_DLC_JITTER64]) < sizeof(jitter64))
                return -1;
            jitter64 = rta_getattr_u64(tb[TCA_DLC_JITTER64]);

        }

        if (tb[TCA_DLC_MODEL]) {
            if (RTA_PAYLOAD(tb[TCA_DLC_MODEL]) < sizeof(*model))
                return -1;
            model = RTA_DATA(tb[TCA_DLC_MODEL]);
        }
        if (tb[TCA_MARKOV_CHAIN]) {
            struct rtattr *lb[MC_STATE_MAX + 1];
            parse_rtattr_nested(lb, MC_STATE_MAX, tb[TCA_MARKOV_CHAIN]);
            if (lb[MC_STATE_SIMPLE]){
                if (RTA_PAYLOAD(tb[MC_STATE_SIMPLE]) < sizeof(*simple_state))
                    return -1;
                simple_state = RTA_DATA(lb[MC_STATE_SIMPLE]);
            }
            if (lb[MC_STATE_QUEUE]){
                if (RTA_PAYLOAD(tb[MC_STATE_QUEUE]) < sizeof(*queue_state))
                    return -1;
                queue_state = RTA_DATA(lb[MC_STATE_QUEUE]);
            }
            if (lb[MC_STATE_LOSS]){
                if (RTA_PAYLOAD(tb[MC_STATE_LOSS]) < sizeof(*loss_state))
                    return -1;
                loss_state = RTA_DATA(lb[MC_STATE_LOSS]);
            }
        }

    }

    open_json_object("dlc_model");
        print_float(PRINT_JSON, "delay", NULL, (double)latency64 / 1000000000.);
        print_float(PRINT_JSON, "jitter", NULL, (double)jitter64 / 1000000000.);
        print_uint(PRINT_JSON, "m/m/1/k rho", NULL, qopt.mm1_rho);
        print_uint(PRINT_JSON, "jitter_steps", NULL, qopt.jitter_steps);
        PRINT_PERCENT("loss", qopt.loss);
        PRINT_PERCENT("mu", qopt.mu);
        print_uint(PRINT_JSON, "mean_burst_len", NULL, qopt.mean_burst_len);
        print_uint(PRINT_JSON, "mean_burst_len", NULL, qopt.mean_good_burst_len);

        rate64 = rate64 ? : qopt.rate;
        tc_print_rate(PRINT_JSON, "rate", " rate %s", rate64);
        print_uint(PRINT_JSON, "limit", NULL, qopt.limit);
    close_json_object();

    if (model == NULL){
        printf("Error: failed to get a dlc_model for print");
        return 0;
    }
    open_json_object("main_chain");
        print_uint(PRINT_JSON, "num_states", NULL, model->mc_num_states);
        print_uint(PRINT_JSON, "curr_state", NULL, model->mc_curr_state);
        print_uint(PRINT_JSON, "delaydist_size", NULL, model->delaydist_size);
        open_json_array(PRINT_JSON, "states");
        if (simple_state != NULL){
            open_json_object("simple_state");
                print_float(PRINT_JSON, "delay", NULL, (double) simple_state->delay / 1000000000.);
                print_float(PRINT_JSON, "jitter", NULL, (double) simple_state->jitter / 1000000000.);
                print_uint(PRINT_JSON, "delaydist_size", NULL, simple_state->delaydist_size);
                print_transition_probs(simple_state->trans_probs);
            close_json_object();
        }
        if (queue_state != NULL){
            open_json_object("queue_state");
                print_uint(PRINT_JSON, "mm1k_num_states", NULL, queue_state->mm1k_num_states);
                print_transition_probs(queue_state->trans_probs);
            close_json_object();
        }
        if (loss_state != NULL){
            open_json_object("loss_state");
                print_float(PRINT_JSON, "delay", NULL, (double) loss_state->max_delay / 1000000000.);
                print_transition_probs(loss_state->trans_probs);
            close_json_object();
        }
        close_json_array(PRINT_JSON, ",\n");
    close_json_object();

    return 0;
}

struct qdisc_util dlc_qdisc_util = {
    .id            = "dlc",
    .parse_qopt    = dlc_parse_opt,
    .print_qopt    = dlc_print_opt,
};
