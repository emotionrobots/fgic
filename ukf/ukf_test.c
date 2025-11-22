#include <stdio.h>
#include <math.h>
#include "ukf.h"

/* Simple deterministic pseudo-random generator for repeatable noise */
static unsigned int lcg_state = 1;
static double lcg_rand_uniform(void)
{
    lcg_state = 1664525u * lcg_state + 1013904223u;
    return (double)(lcg_state & 0xFFFFFF) / (double)0xFFFFFF; /* [0,1) */
}

/* Generate approximate zero-mean noise with magnitude up to std */
static double noise_uniform(double std)
{
    /* Transform uniform [0,1) to [-0.5,0.5) and scale */
    return (lcg_rand_uniform() - 0.5) * 2.0 * std;
}

/* Process model: constant velocity in 1D
 * State: [x, v]^T
 */
static void cv_process_model(double *x, const double *u, double dt)
{
    (void)u; /* unused */
    double pos = x[0];
    double vel = x[1];

    pos = pos + vel * dt;
    /* vel stays constant */

    x[0] = pos;
    x[1] = vel;
}

/* Measurement model: position only */
static void cv_measurement_model(const double *x, double *z)
{
    z[0] = x[0]; /* measure position */
}

static void print_vec(const char *name, const double *v, int n)
{
    printf("%s = [", name);
    for (int i = 0; i < n; ++i) {
        printf(" % .5f", v[i]);
    }
    printf(" ]\n");
}

int main(void)
{
    /* Example UKF for 1D constant velocity */
    const int n_x = 2;
    const int n_z = 1;

    ukf_t ukf;
    ukf_status_t st;

    /* UKF parameters: fairly standard */
    double alpha = 1e-3;
    double beta  = 2.0;
    double kappa = 0.0;

    st = ukf_init(&ukf, n_x, n_z, alpha, beta, kappa);
    if (st != UKF_OK) {
        printf("ukf_init failed: %d\n", st);
        return -1;
    }

    ukf_set_models(&ukf, cv_process_model, cv_measurement_model);

    /* Initial guess: wrong on purpose */
    double x0[2] = {0.0, 0.0};     /* pos=0, vel=0 */
    double P0[4] = {
        10.0, 0.0,
        0.0, 10.0
    }; /* large uncertainty */

    st = ukf_set_state(&ukf, x0, P0);
    if (st != UKF_OK) {
        printf("ukf_set_state failed: %d\n", st);
        return -1;
    }

    /* Process noise: small */
    double q_pos = 1e-4;
    double q_vel = 1e-4;
    double Q[4] = {
        q_pos, 0.0,
        0.0,   q_vel
    };

    /* Measurement noise (normal case) */
    double r_meas = 0.1;
    double R[1] = { r_meas * r_meas };

    st = ukf_set_noise(&ukf, Q, R);
    if (st != UKF_OK) {
        printf("ukf_set_noise failed: %d\n", st);
        return -1;
    }

    /* True state */
    double x_true[2] = {0.0, 1.0};  /* pos=0, vel=1 m/s */

    double dt = 0.1;
    int N = 60;

    printf("Time step dt = %.3f, N = %d\n", dt, N);
    printf("Case 1: Normal operation with measurement every step.\n");
    printf("k, t, x_true, v_true, z_meas, x_est, v_est\n");

    for (int k = 0; k < N; ++k) {
        double t = (k+1) * dt;

        /* Evolve true system (no process noise here, for clarity) */
        x_true[0] = x_true[0] + x_true[1] * dt; /* pos */
        /* vel stays constant */

        /* Generate noisy position measurement */
        double z_meas[1];
        z_meas[0] = x_true[0] + noise_uniform(0.1); /* std dev ~0.1 */

        /* UKF predict-update */
        st = ukf_predict(&ukf, NULL, dt);
        if (st != UKF_OK) {
            printf("ukf_predict failed at k=%d: %d\n", k, st);
            return -1;
        }

        st = ukf_update(&ukf, z_meas);
        if (st != UKF_OK) {
            printf("ukf_update failed at k=%d: %d\n", k, st);
            return -1;
        }

        printf("%2d, %.2f, % .4f, % .4f, % .4f, % .4f, % .4f\n",
               k, t,
               x_true[0], x_true[1],
               z_meas[0],
               ukf.x[0], ukf.x[1]);
    }

    /* Corner case 1: Missing measurements for some steps.
     * We just call predict() without update().
     */

    printf("\nCase 2: Missing measurement from k=20 to k=29.\n");
    /* Reinitialize filter for clarity */
    ukf_set_state(&ukf, x0, P0);
    ukf_set_noise(&ukf, Q, R);

    x_true[0] = 0.0;
    x_true[1] = 1.0;

    for (int k = 0; k < N; ++k) {
        double t = (k+1) * dt;
        x_true[0] = x_true[0] + x_true[1] * dt;

        st = ukf_predict(&ukf, NULL, dt);
        if (st != UKF_OK) {
            printf("ukf_predict failed at k=%d: %d\n", k, st);
            return -1;
        }

        if (k < 20 || k >= 30) {
            double z_meas[1];
            z_meas[0] = x_true[0] + noise_uniform(0.1);
            st = ukf_update(&ukf, z_meas);
            if (st != UKF_OK) {
                printf("ukf_update failed at k=%d: %d\n", k, st);
                return -1;
            }
            printf("k=%2d t=%.2f (with meas) x_true=%.3f x_est=%.3f\n",
                   k, t, x_true[0], ukf.x[0]);
        } else {
            /* No measurement */
            printf("k=%2d t=%.2f (NO meas)  x_true=%.3f x_est=%.3f\n",
                   k, t, x_true[0], ukf.x[0]);
        }
    }

    /* Corner case 2: Very small measurement noise (R ~ 0)
     * Filter will strongly trust measurement and hug z_meas.
     */

    printf("\nCase 3: Very small measurement noise (R ~ 0).\n");
    double R_small[1] = { 1e-8 }; /* nearly noise-free measurement */
    ukf_set_state(&ukf, x0, P0);
    ukf_set_noise(&ukf, Q, R_small);

    x_true[0] = 0.0;
    x_true[1] = 1.0;

    for (int k = 0; k < 10; ++k) {
        double t = (k+1) * dt;
        x_true[0] = x_true[0] + x_true[1] * dt;

        double z_meas[1];
        z_meas[0] = x_true[0] + noise_uniform(0.01); /* tiny noise */

        st = ukf_predict(&ukf, NULL, dt);
        if (st != UKF_OK) {
            printf("ukf_predict failed: %d\n", st);
            return -1;
        }

        st = ukf_update(&ukf, z_meas);
        if (st != UKF_OK) {
            printf("ukf_update failed: %d\n", st);
            return -1;
        }

        printf("k=%2d t=%.2f x_true=%.5f z_meas=%.5f x_est=%.5f\n",
               k, t, x_true[0], z_meas[0], ukf.x[0]);
    }

    return 0;
}

