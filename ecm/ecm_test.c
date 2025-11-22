#include <stdio.h>
#include <math.h>
#include "ecm.h"

/* Simple normal noise generator for measurements */
static unsigned int lcg_state = 123456789u;

static double rand_uniform(void)
{
    lcg_state = 1664525u * lcg_state + 1013904223u;
    return (double)(lcg_state & 0xFFFFFFu) / (double)0xFFFFFFu;
}

static double rand_normal(double std)
{
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    if (u1 < 1e-12) u1 = 1e-12;
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    return std * r * cos(theta);
}

static void print_header(void)
{
    printf("phase,k,t,I,"
           "soc_true,soc_model,"
           "V_true,V_model,"
           "R0_true,R0_model,"
           "R1_true,R1_model,"
           "C1_true,C1_model\n");
}

int main(void)
{
    ecm_t e_true;
    ecm_t e_model;

    ecm_init_default(&e_true);
    ecm_init_default(&e_model);

    /* Make true parameters different from model so adaptation has
     * something to converge to.
     */
    for (int i = 0; i < ECM_TABLE_SIZE; ++i) {
        e_true.r0_table[i] *= 1.5;  /* true cell more resistive */
        e_true.r1_table[i] *= 0.7;  /* faster RC dynamics */
        e_true.c1_table[i] *= 1.4;  /* larger diffusion capacitance */
    }

    /* Reset states */
    ecm_reset_state(&e_true, 0.8, 20.0);
    ecm_reset_state(&e_model, 0.7, 20.0);

    double dt = 1.0;    /* [s] */
    double T_amb = 20.0;
    int k = 0;

    print_header();

    /* Helper macro to run a phase with specified current */
    auto run_phase = [&](int phase_id, double I, int steps) {
        for (int n = 0; n < steps; ++n, ++k) {
            double t = (k + 1) * dt;

            /* True dynamics */
            ecm_step(&e_true, I, T_amb, dt);
            double V_true_clean = ecm_terminal_voltage(&e_true, I);
            double V_meas = V_true_clean + rand_normal(0.005); /* 5 mV noise */

            /* Model dynamics */
            ecm_step(&e_model, I, T_amb, dt);
            /* Use true v_rc as if coming from a perfect UKF */
            double vrc_est = e_true.v_rc;
            ecm_update_from_measurement(&e_model, I, V_meas, vrc_est, dt);

            double V_model = ecm_terminal_voltage(&e_model, I);

            /* Look up effective parameters at current SOC / T */
            double R0_true  = ecm_lookup_r0(&e_true,  e_true.soc,  e_true.T);
            double R0_model = ecm_lookup_r0(&e_model, e_model.soc, e_model.T);
            double R1_true  = ecm_lookup_r1(&e_true,  e_true.soc,  e_true.T);
            double R1_model = ecm_lookup_r1(&e_model, e_model.soc, e_model.T);
            double C1_true  = ecm_lookup_c1(&e_true,  e_true.soc,  e_true.T);
            double C1_model = ecm_lookup_c1(&e_model, e_model.soc, e_model.T);

            printf("%d,%d,%.1f,%.4f,"
                   "%.4f,%.4f,"
                   "%.4f,%.4f,"
                   "%.6f,%.6f,"
                   "%.6f,%.6f,"
                   "%.2f,%.2f\n",
                   phase_id, k, t, I,
                   e_true.soc, e_model.soc,
                   V_true_clean, V_model,
                   R0_true, R0_model,
                   R1_true, R1_model,
                   C1_true, C1_model);
        }
    };

    /* Phase 0: charging at -1 A (raise SOC) */
    run_phase(0, -1.0, 100);

    /* Phase 1: rest after charge */
    run_phase(1, 0.0, 100);

    /* Phase 2: discharging at +1 A */
    run_phase(2, 1.0, 100);

    /* Phase 3: rest after discharge */
    run_phase(3, 0.0, 100);

    return 0;
}

