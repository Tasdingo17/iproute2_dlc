#ifndef _DLC_TCA_SPEC_H
#define _DLC_TCA_SPEC_H

/* TODO: mb latency64, jitter64, mm1_rho64 */
enum {
    TCA_DLC_UNSPEC,
    TCA_DLC_DELAY_DIST,
    TCA_DLC_PAD,
    TCA_DLC_LATENCY64,
    TCA_DLC_JITTER64,
    TCA_DLC_RATE64,
    TCA_DLC_MODEL,      // for dumping
    TCA_MARKOV_CHAIN,   // for dumping
    TCA_MARKOV_PROBS,    // for dumping
    __TCA_DLC_MAX,
};

#define TCA_DLC_MAX (__TCA_DLC_MAX - 1)

struct tc_dlc_qopt {
    __u32   latency;                /* added delay (us) */
    __u32   jitter;                 /* random jitter in latency (us) */
    __u32   mm1_rho;                /* m/m/1/k rho */
    __u32   jitter_steps;           /* number of queue steps (states in m/m/1/k) */
    __u32   loss;                   /* random packet loss (0=none ~0=100%) */
    __u32   mu;                     /* probability to be in queue state */
    __u32   mean_burst_len;         /* mean number of consequitive packet losses */
    __u32   mean_good_burst_len;    /* mean number of consequitive queued packets */
    __u32   limit;                  /* fifo limit (packets) */
    __u32   rate;                   /* max packet send rate */
};

// for dumping
struct tc_dlc_model {
    __u32 mc_num_states;
    __u32 mc_curr_state;
    __u32 delaydist_size;
};

enum {
    MC_STATE_CONST,
    MC_STATE_SIMPLE,
    MC_STATE_QUEUE,
    MC_STATE_LOSS,
    __MC_STATE_MAX
};
#define MC_STATE_MAX (__MC_STATE_MAX - 1)

struct tc_dlc_simple_state{
    __u32 delay;
    __u32 jitter;
    __u32 delaydist_size;
    __u16 trans_probs[3];   /* for specific dlc structure */
};

struct tc_dlc_queue_state{
    __u32 mm1k_num_states;
    __u16 trans_probs[3];   /* for specific dlc structure */
};

struct tc_dlc_loss_state{
    __u32 max_delay;
    __u16 trans_probs[3];   /* for specific dlc structure */
};

#define DLC_DIST_SCALE  8192
#define DLC_DIST_MAX    16384


#endif
