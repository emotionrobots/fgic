/*!
 *======================================================================================================================
 *
 * @file		ukf.h
 *
 * @brief		Unscented Kalman Filter header file
 *
 *======================================================================================================================
 */
#ifndef __UKF_H__
#define __UKF_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Adjust as needed */
#define UKF_MAX_N 		8   /* max state dimension */
#define UKF_MAX_M 		8   /* max measurement dimension */

#define UKF_MAX_SIGMA 		(2 * UKF_MAX_N + 1)


/* Use double for numerical robustness */
typedef double ukf_float;


/* Status codes */
typedef enum 
{
    UKF_OK = 0,
    UKF_ERR_DIM = -1,
    UKF_ERR_CHOL = -2,
    UKF_ERR_INV  = -3
} 
ukf_status_t;


/* Forward declaration */
struct ukf_s;

/*!
 * Process model callback:
 *   x: in/out state vector of size n_x
 *   u: input pointer (may be NULL)
 *   dt: time step
 */
typedef void (*ukf_fx_t)(ukf_float *x, const ukf_float *u, ukf_float dt);


/*! 
 * Measurement model callback:
 *   x: state vector of size n_x
 *   z: output measurement vector of size n_z
 */
typedef void (*ukf_hx_t)(const ukf_float *x, ukf_float *z);


/* Main UKF struct */
typedef struct ukf_s {
    int n_x;   /* state dimension */
    int n_z;   /* measurement dimension */
    int n_sigma; /* 2*n_x + 1 */

    /* State mean and covariance */
    ukf_float x[UKF_MAX_N];
    ukf_float P[UKF_MAX_N * UKF_MAX_N];

    /* Process and measurement noise covariances */
    ukf_float Q[UKF_MAX_N * UKF_MAX_N];
    ukf_float R[UKF_MAX_M * UKF_MAX_M];

    /* Sigma points in state space (predicted) */
    ukf_float sigma_x[UKF_MAX_N * UKF_MAX_SIGMA];

    /* Weights for mean and covariance */
    ukf_float Wm[UKF_MAX_SIGMA];
    ukf_float Wc[UKF_MAX_SIGMA];

    /* Scaling parameters */
    ukf_float alpha;
    ukf_float beta;
    ukf_float kappa;
    ukf_float lambda;

    /* User-supplied models */
    ukf_fx_t fx;
    ukf_hx_t hx;
} 
ukf_t;


ukf_status_t ukf_init(ukf_t *ukf, int n_x, int n_z,
                      ukf_float alpha, ukf_float beta, ukf_float kappa);

ukf_status_t ukf_set_state(ukf_t *ukf,
                           const ukf_float *x0,
                           const ukf_float *P0);

ukf_status_t ukf_set_noise(ukf_t *ukf,
                           const ukf_float *Q,
                           const ukf_float *R);

void ukf_set_models(ukf_t *ukf, ukf_fx_t fx, ukf_hx_t hx);

ukf_status_t ukf_predict(ukf_t *ukf, const ukf_float *u, ukf_float dt);

ukf_status_t ukf_update(ukf_t *ukf, const ukf_float *z);

#ifdef __cplusplus
}
#endif

#endif /* __UKF_H__ */

