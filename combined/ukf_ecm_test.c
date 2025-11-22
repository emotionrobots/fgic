/*
 * ukf_ecm_test.c
 *
 * Combined UKF + ECM unit test:
 *   - "True" ECM cell: e_true (plant)
 *   - "Model" ECM parameters: e_model (used inside UKF process / measurement models)
 *
 * UKF state:
 *   x = [ SOC, VRC, H, T ]^T
 *
 * Measurement:
 *   z = V_term
 *
 * Input to process model:
 *   u = [ I, T_amb ]^T
 *
 * Scenario:
 *   1) Discharging
 *   2) Rest
 *   3) Charging
 *   4) Rest
 *
 * Output: CSV of
 *   k, t, I,
 *   SOC_true, SOC_est,
 *   T_true, T_est,
 *   H_est,
 *   V_true, V_meas, V_est,
 *   R0_true, R0_model,
 *   R1_true, R1_model,
 *   C1_true, C1_model
 *
 * This demonstrates that:
 *   - The UKF tracks SOC and T reasonably
 *   - The ECM + UKF measurement model reproduces terminal voltage
 *   - ECM parameter lookups and Arrhenius scaling are working
 */

#include <stdio.h>
#include <math.h>
#include "ukf.h"
#include "ecm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* ---------------- Random noise generator (for V_meas) ----------------- */

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

/* ---------------- Globals: ECM model used by UKF, current for hx ------- */

static ecm_t g_ecm_model;     /* parameter container for UKF model */
static double g_current = 0.0; /* current used in measurement model */

/* ---------------- UKF Process Model: fx(x,u,dt) ------------------------
 *
 * State: x = [ SOC, VRC, H, T ]
 * Input: u = [ I, T_amb ]
 *
 * Uses g_ecm_model tables (R0/R1/C1, H_chg/H_dsg) via ecm_* lookups.
 * Does NOT modify g_ecm_model dynamic state (we treat it as parameter store).
 */

static void ukf_ecm_process_model(double *x, const double *u, double dt)
{
    double soc = x[0];
    double vrc = x[1];
    double H   = x[2];
    double T   = x[3];

    double I     = u[0];
    double T_amb = u[1];

    /* Capacity in Coulombs from model ECM (Q_Ah) */
    double Qc = g_ecm_model.Q_Ah * 3600.0;

    /* SOC dynamics: I>0 = discharge */
    double dSOC = -(I * dt) / Qc;
    soc += dSOC;
    if (soc < 0.0) soc = 0.0;
    if (soc > 1.0) soc = 1.0;

    /* RC branch */
    double R1 = ecm_lookup_r1(&g_ecm_model, soc, T);
    double C1 = ecm_lookup_c1(&g_ecm_model, soc, T);
    double tau1 = R1 * C1;
    if (tau1 < 1e-9) tau1 = 1e-9;

    double dVRC = dt * ( -vrc / tau1 + I / C1 );
    vrc += dVRC;

    /* Hysteresis dynamics: first-order towards table hysteresis */
    const double I_eps = 1e-3;
    int dir = 0;
    if (I > I_eps)      dir = +1; /* discharge */
    else if (I < -I_eps) dir = -1; /* charge */

    double H_inf = 0.0;
    if (dir > 0) {
        H_inf = ecm_lookup_h_dsg(&g_ecm_model, soc);
    } else if (dir < 0) {
        H_inf = ecm_lookup_h_chg(&g_ecm_model, soc);
    } else {
        /* near rest: keep H roughly where it is (small relaxation) */
        H_inf = H;
    }

    double tau_H = 200.0; /* example hysteresis time constant [s] */
    double dH = dt * (H_inf - H) / tau_H;
    H += dH;

    /* Thermal dynamics */
    double R0 = ecm_lookup_r0(&g_ecm_model, soc, T);
    double power_loss = I * I * R0;
    double dT = dt * ( (power_loss - (T - T_amb) / g_ecm_model.R_th) / g_ecm_model.C_th );
    T += dT;

    /* Clamp T to reasonable bounds */
    if (T < -40.0) T = -40.0;
    if (T > 80.0)  T = 80.0;

    x[0] = soc;
    x[1] = vrc;
    x[2] = H;
    x[3] = T;
}

/* ---------------- UKF Measurement Model: hx(x) -------------------------
 *
 * z = [ V_term ]
 *
 * Uses g_ecm_model lookups and global g_current.
 */

static void ukf_ecm_measurement_model(const double *x, double *z)
{
    double soc = x[0];
    double vrc = x[1];
    double H   = x[2];
    double T   = x[3];

    double I = g_current;

    double ocv = ecm_lookup_ocv(&g_ecm_model, soc, T);
    double R0  = ecm_lookup_r0(&g_ecm_model, soc, T);

    double V = ocv + H + vrc - I * R0;
    z[0] = V;
}

/* ---------------- Helper: CSV header and line -------------------------- */

static void print_header(void)
{
    printf("k,t,I,"
           "SOC_true,SOC_est,"
           "T_true,T_est,"
           "H_est,"
           "V_true,V_meas,V_est,"
           "R0_true,R0_model,"
           "R1_true,R1_model,"
           "C1_true,C1_model\n");
}

static void print_line(int k, double t, double I,
                       const ecm_t *e_true,
                       const ukf_t *ukf,
                       double V_true_clean,
                       double V_meas,
                       double V_est)
{
    double soc_true = e_true->soc;
    double T_true   = e_true->T;

    double soc_est  = ukf->x[0];
    double vrc_est  = ukf->x[1];
    double H_est    = ukf->x[2];
    double T_est    = ukf->x[3];

    /* Parameter lookups */
    double R0_true  = ecm_lookup_r0(e_true,  soc_true, T_true);
    double R1_true  = ecm_lookup_r1(e_true,  soc_true, T_true);
    double C1_true  = ecm_lookup_c1(e_true,  soc_true, T_true);

    double R0_model = ecm_lookup_r0(&g_ecm_model, soc_est, T_est);
    double R1_model = ecm_lookup_r1(&g_ecm_model, soc_est, T_est);
    double C1_model = ecm_lookup_c1(&g_ecm_model, soc_est, T_est);

    printf("%d,%.2f,%.4f,"
           "%.5f,%.5f,"
           "%.2f,%.2f,"
           "%.5f,"
           "%.4f,%.4f,%.4f,"
           "%.6f,%.6f,"
           "%.6f,%.6f,"
           "%.2f,%.2f\n",
           k, t, I,
           soc_true, soc_est,
           T_true, T_est,
           H_est,
           V_true_clean, V_meas, V_est,
           R0_true, R0_model,
           R1_true, R1_model,
           C1_true, C1_model);
}

/* ---------------- Main test ------------------------------------------- */

int main(void)
{
    /* True ECM and model ECM */
    ecm_t e_true;
    ecm_t e_model;

    ecm_init_default(&e_true);
    ecm_init_default(&e_model);

    /* Make the true cell different so UKF+ECM have something to learn */
    for (int i = 0; i < ECM_TABLE_SIZE; ++i) {
        e_true.r0_table[i] *= 1.4;  /* more resistive */
        e_true.r1_table[i] *= 0.8;  /* slightly smaller R1 */
        e_true.c1_table[i] *= 1.3;  /* larger diffusion cap */
    }

    /* Initial states */
    ecm_reset_state(&e_true, 0.8, 20.0);  /* true cell at SOC=0.8, T=20C */
    ecm_reset_state(&e_model, 0.7, 20.0); /* model parameters (not used as dynamics) */

    /* Copy model ECM into global for UKF */
    g_ecm_model = e_model;

    /* UKF setup */
    const int n_x = 4; /* SOC, VRC, H, T */
    const int n_z = 1; /* terminal voltage */

    ukf_t ukf;
    ukf_status_t st;

    double alpha = 1e-3;
    double beta  = 2.0;
    double kappa = 0.0;

    st = ukf_init(&ukf, n_x, n_z, alpha, beta, kappa);
    if (st != UKF_OK) {
        printf("ukf_init failed: %d\n", st);
        return -1;
    }

    ukf_set_models(&ukf, ukf_ecm_process_model, ukf_ecm_measurement_model);

    /* Initial UKF state: deliberately a bit wrong */
    double x0[4] = {
        0.6,   /* SOC */
        0.0,   /* VRC */
        0.0,   /* H */
        25.0   /* T (a bit off) */
    };

    /* Initial covariance */
    double P0[16] = {
        0.04, 0,    0,    0,
        0,    0.01, 0,    0,
        0,    0,    0.01, 0,
        0,    0,    0,    4.0
    };

    st = ukf_set_state(&ukf, x0, P0);
    if (st != UKF_OK) {
        printf("ukf_set_state failed: %d\n", st);
        return -1;
    }

    /* Process noise Q */
    double Q[16] = {0};
    Q[0*4 + 0] = 1e-6; /* SOC */
    Q[1*4 + 1] = 1e-5; /* VRC */
    Q[2*4 + 2] = 1e-6; /* H */
    Q[3*4 + 3] = 5e-4; /* T */

    /* Measurement noise R (10 mV) */
    double R[1] = { 0.01 * 0.01 };

    st = ukf_set_noise(&ukf, Q, R);
    if (st != UKF_OK) {
        printf("ukf_set_noise failed: %d\n", st);
        return -1;
    }

    /* Simulation setup */
    double dt    = 1.0;
    double T_amb = 20.0;
    int    N     = 400;

    print_header();

    for (int k = 0; k < N; ++k) {
        double t = (k + 1) * dt;

        /* Current profile:
         *  0-149 : discharge +1.0 A
         *  150-199: rest (0 A)
         *  200-299: charge -0.5 A
         *  300-399: rest (0 A)
         */
        double I;
        if (k < 150)        I =  1.0;
        else if (k < 200)   I =  0.0;
        else if (k < 300)   I = -0.5;
        else                I =  0.0;

        /* ---- True ECM dynamics ---- */
        ecm_step(&e_true, I, T_amb, dt);
        double V_true_clean = ecm_terminal_voltage(&e_true, I);

        /* Add measurement noise */
        double V_meas = V_true_clean + rand_normal(0.01);

        /* ---- UKF predict ---- */
        double u[2] = { I, T_amb };
        st = ukf_predict(&ukf, u, dt);
        if (st != UKF_OK) {
            printf("ukf_predict failed at k=%d: %d\n", k, st);
            return -1;
        }

        /* ---- UKF update ---- */
        g_current = I; /* set global current for hx */
        double z[1] = { V_meas };
        st = ukf_update(&ukf, z);
        if (st != UKF_OK) {
            printf("ukf_update failed at k=%d: %d\n", k, st);
            return -1;
        }

        /* ---- ECM parameter and hysteresis updates from UKF ---- */
        {
            double soc_est = ukf.x[0];
            double vrc_est = ukf.x[1];
            double H_est   = ukf.x[2];

            /* Feed UKF's VRC into ECM parameter learner */
            ecm_update_from_measurement(&g_ecm_model, I, V_meas, vrc_est, dt);

            /* Update hysteresis tables from UKF H_est */
            int is_chg = (I < -1e-3) ? 1 : 0;
            ecm_update_h_from_ukf(&g_ecm_model, soc_est, H_est, is_chg);
        }

        /* Predicted measurement from UKF model (noise-free) */
        double z_est[1];
        ukf_ecm_measurement_model(ukf.x, z_est);
        double V_est = z_est[0];

        print_line(k, t, I, &e_true, &ukf, V_true_clean, V_meas, V_est);
    }

    return 0;
}

