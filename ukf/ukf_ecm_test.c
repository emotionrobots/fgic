/*
 * ukf_ecm_test.c
 *
 * UKF unit test with ECM-like system model:
 *   State x = [ SOC, VRC, H, T ]^T
 *
 * SOC : State of Charge (0..1)
 * VRC : RC branch voltage (V)
 * H   : Hysteresis voltage (V)
 * T   : Cell core temperature (degC)
 *
 * Input u = [ I, T_amb ]^T:
 *   I     : Cell current (A), I > 0 = discharge
 *   T_amb : Ambient temperature (degC)
 *
 * Measurement:
 *   z = V_term (terminal voltage)
 *
 * The test simulates a "true" battery with some parameters and runs a UKF
 * with slightly mismatched parameters to estimate SOC, VRC, H, T from
 * noisy terminal voltage measurements.
 */

#include <stdio.h>
#include <math.h>
#include "ukf.h"

/* ---------------- True ECM Parameters (simulation model) ---------------- */

static const double Q_Ah_true  = 2.5;     /* Cell capacity [Ah] */
static const double R0_true    = 0.020;   /* Ohmic resistance [ohm] at 25C */
static const double R1_true    = 0.005;   /* RC branch resistance [ohm] */
static const double C1_true    = 500.0;   /* RC branch capacitance [F] */
static const double H0_true    = 0.03;    /* Max hysteresis [V] */
static const double kH_true    = 2.0;     /* Hysteresis nonlinearity */
static const double tauH_true  = 200.0;   /* Hysteresis time constant [s] */
static const double Cth_true   = 200.0;   /* Thermal capacity [J/degC] */
static const double Rth_true   = 3.0;     /* Thermal resistance [degC/W] */
static const double alpha_R0_true = 0.01; /* Temp coefficient for R0 [1/degC] */

/* UKF will use slightly mismatched parameters */
static const double Q_Ah_filt  = 2.3;
static const double R0_filt    = 0.022;
static const double R1_filt    = 0.006;
static const double C1_filt    = 450.0;
static const double H0_filt    = 0.025;
static const double kH_filt    = 2.0;
static const double tauH_filt  = 180.0;
static const double Cth_filt   = 180.0;
static const double Rth_filt   = 3.0;
static const double alpha_R0_filt = 0.01;

/* Simple OCV vs SOC function (same for true and filter here).
 * For a real system you'd replace this with a lookup + interpolation.
 */
static double ocv_from_soc(double soc)
{
    /* Clamp SOC */
    if (soc < 0.0) soc = 0.0;
    if (soc > 1.0) soc = 1.0;

    /* Very simple model: 3.0 V at 0% SOC, 3.7 V at 100% SOC */
    return 3.0 + 0.7 * soc;
}

/* Pseudo-random generator (LCG) to create repeatable noise */
static unsigned int lcg_state = 1234567u;

static double rand_uniform(void)
{
    lcg_state = 1664525u * lcg_state + 1013904223u;
    return (double)(lcg_state & 0xFFFFFFu) / (double)0xFFFFFFu;
}

static double rand_normal(double std)
{
    /* Box-Muller transform */
    double u1 = rand_uniform();
    double u2 = rand_uniform();

    if (u1 < 1e-12) u1 = 1e-12; /* avoid log(0) */

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    double n = r * cos(theta);  /* standard normal */
    return std * n;
}

/* ---------------- ECM truth model (for simulation) ---------------------- */

/* Simulate one step of the "true" ECM model.
 * state_true = [ SOC, VRC, H, T ].
 * u = [ I, T_amb ].
 */
static void ecm_true_step(double *x_true, const double *u, double dt)
{
    double soc = x_true[0];
    double vrc = x_true[1];
    double H   = x_true[2];
    double T   = x_true[3];

    double I     = u[0];
    double T_amb = u[1];

    /* Capacity in Coulombs */
    double Qc = Q_Ah_true * 3600.0;

    /* Temperature-dependent R0 (linear model) */
    double R0_T = R0_true * (1.0 + alpha_R0_true * (T - 25.0));

    /* SOC dynamics (I > 0 is discharge) */
    double dSOC = - (I * dt) / Qc;
    soc += dSOC;

    /* RC branch */
    double tau1 = R1_true * C1_true; /* time constant */
    double dVRC = dt * ( -vrc / tau1 + I / C1_true );
    vrc += dVRC;

    /* Hysteresis (first-order towards H_inf) */
    double H_inf = H0_true * tanh(kH_true * I);
    double dH = dt * ( (H_inf - H) / tauH_true );
    H += dH;

    /* Thermal dynamics */
    double power_loss = I * I * R0_T; /* Joule heating */
    double dT = dt * ( (power_loss - (T - T_amb) / Rth_true) / Cth_true );
    T += dT;

    /* Simple clamping */
    if (soc < 0.0) soc = 0.0;
    if (soc > 1.0) soc = 1.0;
    if (T < -40.0) T = -40.0;
    if (T > 80.0)  T = 80.0;

    x_true[0] = soc;
    x_true[1] = vrc;
    x_true[2] = H;
    x_true[3] = T;
}

/* Compute true terminal voltage given state_true and input u = [I, T_amb]. */
static double ecm_true_voltage(const double *x_true, const double *u)
{
    double soc = x_true[0];
    double vrc = x_true[1];
    double H   = x_true[2];
    double T   = x_true[3];

    double I = u[0];

    double R0_T = R0_true * (1.0 + alpha_R0_true * (T - 25.0));
    double ocv = ocv_from_soc(soc);

    double V = ocv + H + vrc - I * R0_T;
    return V;
}

/* ---------------- UKF Model (fx, hx) ----------------------------------- */

/* Global current used in measurement model (UKF hx).
 * We update this from main() before each ukf_update().
 */
static double g_meas_current = 0.0;

/* Process model for the UKF.
 * x = [ SOC, VRC, H, T ]
 * u = [ I, T_amb ]
 */
static void ecm_ukf_process_model(double *x, const double *u, double dt)
{
    double soc = x[0];
    double vrc = x[1];
    double H   = x[2];
    double T   = x[3];

    double I     = u[0];
    double T_amb = u[1];

    double Qc = Q_Ah_filt * 3600.0;
    double R0_T = R0_filt * (1.0 + alpha_R0_filt * (T - 25.0));

    /* SOC update */
    double dSOC = - (I * dt) / Qc;
    soc += dSOC;

    /* RC branch */
    double tau1 = R1_filt * C1_filt;
    double dVRC = dt * ( -vrc / tau1 + I / C1_filt );
    vrc += dVRC;

    /* Hysteresis */
    double H_inf = H0_filt * tanh(kH_filt * I);
    double dH = dt * ( (H_inf - H) / tauH_filt );
    H += dH;

    /* Thermal */
    double power_loss = I * I * R0_T;
    double dT = dt * ( (power_loss - (T - T_amb) / Rth_filt) / Cth_filt );
    T += dT;

    /* Clamp */
    if (soc < 0.0) soc = 0.0;
    if (soc > 1.0) soc = 1.0;
    if (T < -40.0) T = -40.0;
    if (T > 80.0)  T = 80.0;

    x[0] = soc;
    x[1] = vrc;
    x[2] = H;
    x[3] = T;
}

/* Measurement model for the UKF.
 * x = [ SOC, VRC, H, T ]
 * Measurement: z = terminal voltage, using same form as truth but with
 * filter parameters and the *known* current g_meas_current.
 */
static void ecm_ukf_measurement_model(const double *x, double *z)
{
    double soc = x[0];
    double vrc = x[1];
    double H   = x[2];
    double T   = x[3];

    double I = g_meas_current;

    double R0_T = R0_filt * (1.0 + alpha_R0_filt * (T - 25.0));
    double ocv = ocv_from_soc(soc);

    double V = ocv + H + vrc - I * R0_T;

    z[0] = V;
}

/* ---------------- Utility Printing ------------------------------------- */

static void print_header(void)
{
    printf("k, t_s, I_A, SOC_true, SOC_est, V_true, V_meas, V_est, T_true, T_est\n");
}

static void print_row(int k, double t, double I,
                      const double *x_true, const ukf_t *ukf,
                      double V_true, double V_meas, double V_est)
{
    printf("%4d, %7.2f, %6.3f, %7.4f, %7.4f, %7.4f, %7.4f, %7.4f, %7.3f, %7.3f\n",
           k, t, I,
           x_true[0], ukf->x[0],
           V_true, V_meas, V_est,
           x_true[3], ukf->x[3]);
}

/* ---------------- Main Unit Test --------------------------------------- */

int main(void)
{
    const int n_x = 4; /* SOC, VRC, H, T */
    const int n_z = 1; /* terminal voltage */

    ukf_t ukf;
    ukf_status_t st;

    /* UKF parameters */
    double alpha = 1e-3;
    double beta  = 2.0;
    double kappa = 0.0;

    st = ukf_init(&ukf, n_x, n_z, alpha, beta, kappa);
    if (st != UKF_OK) {
        printf("ukf_init failed: %d\n", st);
        return -1;
    }

    ukf_set_models(&ukf, ecm_ukf_process_model, ecm_ukf_measurement_model);

    /* Initial UKF state (deliberately off) */
    double x0[4] = {
        0.8,   /* SOC guess */
        0.0,   /* VRC guess */
        0.0,   /* H guess */
        25.0   /* T guess */
    };

    /* Initial covariance: fairly large uncertainty */
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

    /* Process noise covariance Q (diagonal) */
    double Q[16] = {0};
    Q[0*4 + 0] = 1e-6;  /* SOC */
    Q[1*4 + 1] = 1e-5;  /* VRC */
    Q[2*4 + 2] = 1e-6;  /* H */
    Q[3*4 + 3] = 1e-3;  /* T */

    /* Measurement noise covariance R */
    double R[1] = { (0.01 * 0.01) }; /* 10 mV std dev */

    st = ukf_set_noise(&ukf, Q, R);
    if (st != UKF_OK) {
        printf("ukf_set_noise failed: %d\n", st);
        return -1;
    }

    /* True state initialization */
    double x_true[4] = {
        0.9,   /* true SOC */
        0.0,   /* true VRC */
        0.0,   /* true hysteresis */
        25.0   /* true T */
    };

    /* Simulation parameters */
    double dt = 1.0;     /* 1 second time step */
    int N = 300;         /* total steps */
    double T_amb = 25.0; /* constant ambient */

    print_header();

    for (int k = 0; k < N; ++k) {
        double t = (k + 1) * dt;

        /* Simple current profile:
         *  0-99:   1.0 A discharge
         *  100-149: rest
         *  150-199: -0.5 A charge
         *  200-249: 1.5 A discharge
         *  250-299: rest
         */
        double I;
        if (k < 100)        I = 1.0;
        else if (k < 150)   I = 0.0;
        else if (k < 200)   I = -0.5;
        else if (k < 250)   I = 1.5;
        else                I = 0.0;

        double u[2] = { I, T_amb };

        /* --- Simulate true ECM --- */
        ecm_true_step(x_true, u, dt);
        double V_true = ecm_true_voltage(x_true, u);

        /* Add measurement noise */
        double V_meas = V_true + rand_normal(0.01);

        /* --- UKF predict --- */
        st = ukf_predict(&ukf, u, dt);
        if (st != UKF_OK) {
            printf("ukf_predict failed at k=%d: %d\n", k, st);
            return -1;
        }

        /* Update global current for measurement model */
        g_meas_current = I;

        /* --- UKF update --- */
        double z[1] = { V_meas };
        st = ukf_update(&ukf, z);
        if (st != UKF_OK) {
            printf("ukf_update failed at k=%d: %d\n", k, st);
            return -1;
        }

        /* For logging, compute UKF's predicted voltage (noise-free) */
        double V_est;
        {
            double z_est[1];
            ecm_ukf_measurement_model(ukf.x, z_est);
            V_est = z_est[0];
        }

        /* Print every step (you can reduce to every 5 or 10 for less output) */
        print_row(k, t, I, x_true, &ukf, V_true, V_meas, V_est);
    }

    return 0;
}

