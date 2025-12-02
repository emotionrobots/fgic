/*!
 *======================================================================================================================
 *
 * @file	ukf.c
 *
 * @brief	Unscented Kalman Filter implementation
 *
 *======================================================================================================================
 */
#include "ukf.h"
#include <math.h>
#include <string.h>


/* Index helper: row-major */
#define IDX(r,c,nc) 		((r)*(nc) + (c))		


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_zero(ukf_float *A, int rows, int cols)
 *
 *  @brief  	Zero out a matrix
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_zero(ukf_float *A, int rows, int cols)
{
    int n = rows * cols;
    for (int i = 0; i < n; ++i) A[i] = 0.0;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_copy(ukf_float *dst, const ukf_float *src, int rows, int cols)
 *  
 *  @breif	Copy a matrix
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_copy(ukf_float *dst, const ukf_float *src, int rows, int cols)
{
    memcpy(dst, src, sizeof(ukf_float) * rows * cols);
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_add(ukf_float *C, const ukf_float *A, const ukf_float *B, int rows, int cols)
 *
 *  @brief	Add two matrix C = A + B
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_add(ukf_float *C, const ukf_float *A, const ukf_float *B, int rows, int cols)
{
    int n = rows * cols;
    for (int i = 0; i < n; ++i) C[i] = A[i] + B[i];
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_sub(ukf_float *C, const ukf_float *A, const ukf_float *B, int rows, int cols)
 *
 *  @brief	Subtract two matrix:  C = A - B
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_sub(ukf_float *C, const ukf_float *A, const ukf_float *B, int rows, int cols)
{
    int n = rows * cols;
    for (int i = 0; i < n; ++i) C[i] = A[i] - B[i];
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_mul(ukf_float *C, const ukf_float *A, const ukf_float *B, int a_rows, int a_cols, int b_cols)
 *  
 *  @brief	Multiply two matricies: C = A * B
 *
 *  @note	C = A(a_rows x a_cols) * B(a_cols x b_cols) 
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_mul(ukf_float *C, const ukf_float *A, const ukf_float *B, int a_rows, int a_cols, int b_cols)
{
    /* C = A(a_rows x a_cols) * B(a_cols x b_cols) */
    mat_zero(C, a_rows, b_cols);
    for (int i = 0; i < a_rows; ++i) {
        for (int k = 0; k < a_cols; ++k) {
            ukf_float aik = A[IDX(i,k,a_cols)];
            for (int j = 0; j < b_cols; ++j) {
                C[IDX(i,j,b_cols)] += aik * B[IDX(k,j,b_cols)];
            }
        }
    }
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void mat_transpose(ukf_float *AT, const ukf_float *A, int rows, int cols)
 *
 *  @brief	Transpose a matrix  AT = transpose(A)
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void mat_transpose(ukf_float *AT, const ukf_float *A, int rows, int cols)
{
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            AT[IDX(c,r,rows)] = A[IDX(r,c,cols)];
        }
    }
}


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_copy(ukf_float *dst, const ukf_float *src, int n)
 *
 *  @brief	Copy a vector
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void vec_copy(ukf_float *dst, const ukf_float *src, int n)
{
    memcpy(dst, src, sizeof(ukf_float) * n);
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void vec_zero(ukf_float *v, int n)
 *
 *  @brief	Zero a vector
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void vec_zero(ukf_float *v, int n)
{
    for (int i = 0; i < n; ++i) v[i] = 0.0;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 * @fn		void vec_axpy(ukf_float *y, ukf_float a, const ukf_float *x, int n)
 *
 * @brief	Multiply and accumulate
 *
 * @note	y += a * x 
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void vec_axpy(ukf_float *y, ukf_float a, const ukf_float *x, int n)
{
    for (int i = 0; i < n; ++i) 
    {
        y[i] += a * x[i];
    }
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		int chol_decomp(const ukf_float *A, ukf_float *L, int n)
 *
 *  @brief	Cholesky decomposition for SPD matrix A (n x n).
 *
 *  @param	A: input (not modified)
 *  @param	L: output lower-triangular (n x n)
 *
 *  @returns 	0 on success, non-zero on failure (not SPD)
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
int chol_decomp(const ukf_float *A, ukf_float *L, int n)
{
    /* Copy A to L, then factor in place (lower) */
    mat_copy(L, A, n, n);

    for (int i = 0; i < n; ++i) 
    {
        ukf_float sum = L[IDX(i,i,n)];
        for (int k = 0; k < i; ++k) 
	{
            ukf_float Lik = L[IDX(i,k,n)];
            sum -= Lik * Lik;
        }

        if (sum <= 0.0) 
	{
            return -1; /* Not SPD */
        }

        ukf_float diag = sqrt(sum);
        L[IDX(i,i,n)] = diag;

        for (int j = i+1; j < n; ++j) 
	{
            ukf_float s = L[IDX(j,i,n)];
            for (int k = 0; k < i; ++k) 
	    {
                s -= L[IDX(j,k,n)] * L[IDX(i,k,n)];
            }
            L[IDX(j,i,n)] = s / diag;
        }

        /* Zero the upper triangle */
        for (int j = i+1; j < n; ++j) 
	{
            L[IDX(i,j,n)] = 0.0;
        }
    }
    return 0;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		int mat_inv(const ukf_float *A, ukf_float *A_inv, int n)
 *
 *  @brief	Gauss-Jordan matrix inversion with partial pivoting.
 * 
 *  @param	A: input matrix (n x n)
 *  @param	A_inv: output inverse
 * 
 *  @returns 	0 on success, non-zero on failure (singular).
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
int mat_inv(const ukf_float *A, ukf_float *A_inv, int n)
{
    /* Augment [A | I] and perform elimination in place on a temp buffer. */
    ukf_float tmp[UKF_MAX_M * 2 * UKF_MAX_M]; /* fits up to 8x8 */
    if (n > UKF_MAX_M) return -1;

    /* Build augmented matrix */
    for (int r = 0; r < n; ++r) 
    {
        for (int c = 0; c < n; ++c) 
	{
            tmp[IDX(r, c, 2*n)] = A[IDX(r, c, n)];
        }

        for (int c = 0; c < n; ++c) 
	{
            tmp[IDX(r, c+n, 2*n)] = (r == c) ? 1.0 : 0.0;
        }
    }

    /* Gauss-Jordan */
    for (int i = 0; i < n; ++i) 
    {
        /* Pivot search */
        int pivot = i;
        ukf_float max_val = fabs(tmp[IDX(i, i, 2*n)]);
        for (int r = i+1; r < n; ++r) 
	{
            ukf_float v = fabs(tmp[IDX(r, i, 2*n)]);
            if (v > max_val) {
                max_val = v;
                pivot = r;
            }
        }

        if (max_val < 1e-15) 
	{
            return -1; /* Singular */
        }

        /* Swap rows if needed */
        if (pivot != i) 
	{
            for (int c = 0; c < 2*n; ++c) 
	    {
                ukf_float tmpv = tmp[IDX(i,c,2*n)];
                tmp[IDX(i,c,2*n)] = tmp[IDX(pivot,c,2*n)];
                tmp[IDX(pivot,c,2*n)] = tmpv;
            }
        }

        /* Normalize pivot row */
        ukf_float piv = tmp[IDX(i,i,2*n)];
        for (int c = 0; c < 2*n; ++c) 
	{
            tmp[IDX(i,c,2*n)] /= piv;
        }

        /* Eliminate others */
        for (int r = 0; r < n; ++r) 
	{
            if (r == i) continue;
            ukf_float factor = tmp[IDX(r,i,2*n)];
            for (int c = 0; c < 2*n; ++c) 
	    {
                tmp[IDX(r,c,2*n)] -= factor * tmp[IDX(i,c,2*n)];
            }
        }
    }

    /* Extract inverse */
    for (int r = 0; r < n; ++r) 
    {
        for (int c = 0; c < n; ++c) 
	{
            A_inv[IDX(r,c,n)] = tmp[IDX(r,c+n,2*n)];
        }
    }
    return 0;
}


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void ukf_compute_weights(ukf_t *ukf)
 *
 *  @brief	Initialize weights given n_x and alpha, beta, kappa 
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
void ukf_compute_weights(ukf_t *ukf)
{
    int n = ukf->n_x;
    ukf->lambda = ukf->alpha * ukf->alpha * (n + ukf->kappa) - n;
    ukf_float lambda = ukf->lambda;
    ukf_float c = n + lambda;

    ukf->n_sigma = 2 * n + 1;

    ukf->Wm[0] = lambda / c;
    ukf->Wc[0] = lambda / c + (1.0 - ukf->alpha * ukf->alpha + ukf->beta);

    for (int i = 1; i < ukf->n_sigma; ++i) 
    {
        ukf->Wm[i] = 1.0 / (2.0 * c);
        ukf->Wc[i] = 1.0 / (2.0 * c);
    }
}

/* Public API */

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_init(ukf_t *ukf, int n_x, int n_z, ukf_float alpha, ukf_float beta, ukf_float kappa)
 *
 *  @brief 	Initialize UKF dimensions and scaling parameters.
 *
 *  @param	n_x: state dimension (<= UKF_MAX_N)
 *  @param	n_z: measurement dimension (<= UKF_MAX_M)
 *  @param	alpha, beta, kappa: UKF scaling parameters
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
ukf_status_t ukf_init(ukf_t *ukf, 
		      int n_x, int n_z,
                      ukf_float alpha, ukf_float beta, ukf_float kappa)
{
    if (!ukf) return UKF_ERR_DIM;
    if (n_x <= 0 || n_x > UKF_MAX_N) return UKF_ERR_DIM;
    if (n_z <= 0 || n_z > UKF_MAX_M) return UKF_ERR_DIM;

    ukf->n_x = n_x;
    ukf->n_z = n_z;
    ukf->alpha = alpha;
    ukf->beta = beta;
    ukf->kappa = kappa;
    ukf_compute_weights(ukf);

    vec_zero(ukf->x, n_x);
    mat_zero(ukf->P, n_x, n_x);
    mat_zero(ukf->Q, n_x, n_x);
    mat_zero(ukf->R, n_z, n_z);
    mat_zero(ukf->sigma_x, n_x, ukf->n_sigma);

    ukf->fx = NULL;
    ukf->hx = NULL;

    return UKF_OK;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_set_state(ukf_t *ukf, const ukf_float *x0, const ukf_float *P0)
 *
 *  @brief	Set state and covariance 
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
ukf_status_t ukf_set_state(ukf_t *ukf,
                           const ukf_float *x0,
                           const ukf_float *P0)
{
    if (!ukf || !x0 || !P0) return UKF_ERR_DIM;
    int n = ukf->n_x;

    vec_copy(ukf->x, x0, n);
    mat_copy(ukf->P, P0, n, n);

    return UKF_OK;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_set_noise(ukf_t *ukf, const ukf_float *Q, const ukf_float *R)
 *
 *  @brief	Set process and measurement noise covariance 
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
ukf_status_t ukf_set_noise(ukf_t *ukf, const ukf_float *Q, const ukf_float *R)
{
    if (!ukf || !Q || !R) return UKF_ERR_DIM;
    mat_copy(ukf->Q, Q, ukf->n_x, ukf->n_x);
    mat_copy(ukf->R, R, ukf->n_z, ukf->n_z);
    return UKF_OK;
}


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void ukf_set_models(ukf_t *ukf, ukf_fx_t fx, ukf_hx_t hx)
 *
 *  @brief	Set process and measurement model callbacks
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
void ukf_set_models(ukf_t *ukf, ukf_fx_t fx, ukf_hx_t hx)
{
    ukf->fx = fx;
    ukf->hx = hx;
}


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_sigma_points(const ukf_t *ukf, ukf_float *X)
 *
 *  @brief	Generate sigma points from current (x, P) into local X (n_x x n_sigma).
 *              Also scales P by (n_x + lambda) before chol.
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
static 
ukf_status_t ukf_sigma_points(const ukf_t *ukf, ukf_float *X)
{
    int n = ukf->n_x;
    int ns = ukf->n_sigma;
    ukf_float lambda = ukf->lambda;
    ukf_float c = n + lambda;

    /* X[:,0] = x */
    for (int i = 0; i < n; ++i) 
    {
        X[IDX(i,0,ns)] = ukf->x[i];
    }

    /* Compute sqrt(c * P) via Cholesky */
    ukf_float A[UKF_MAX_N * UKF_MAX_N];
    ukf_float L[UKF_MAX_N * UKF_MAX_N];

    /* A = c * P */
    for (int i = 0; i < n*n; ++i) 
    {
        A[i] = c * ukf->P[i];
    }

    if (chol_decomp(A, L, n) != 0) 
    {
        /* Try small jitter for numerical robustness */
        const ukf_float jitter = 1e-9;
        for (int i = 0; i < n; ++i) 
	{
            A[IDX(i,i,n)] += jitter;
        }
        if (chol_decomp(A, L, n) != 0) 
	{
            return UKF_ERR_CHOL;
        }
    }

    /* Remaining sigma points */
    for (int col = 0; col < n; ++col) 
    {
        for (int i = 0; i < n; ++i) 
	{
            ukf_float delta = L[IDX(i,col,n)];
            X[IDX(i, 1 + col, ns)]     = ukf->x[i] + delta;
            X[IDX(i, 1 + col + n, ns)] = ukf->x[i] - delta;
        }
    }

    return UKF_OK;
}

/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_predict(ukf_t *ukf, const ukf_float *u, ukf_float dt)
 *
 *  @brief	Prediction step: propagates state and covariance through fx.
 *  
 *  @param	u: input vector pointer (may be NULL if fx ignores it)
 *  @param	dt: timestep
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
ukf_status_t ukf_predict(ukf_t *ukf, const ukf_float *u, ukf_float dt)
{
    if (!ukf || !ukf->fx) return UKF_ERR_DIM;
    int n = ukf->n_x;
    int ns = ukf->n_sigma;

    ukf_float X[UKF_MAX_N * UKF_MAX_SIGMA];

    /* 1. Generate sigma points from current x, P */
    ukf_status_t st = ukf_sigma_points(ukf, X);
    if (st != UKF_OK) return st;

    /* 2. Propagate through process model: sigma_x = f(X) */
    for (int k = 0; k < ns; ++k) 
    {
        /* Copy column k to temporay state, apply fx in-place, then store */
        ukf_float xk[UKF_MAX_N];
        for (int i = 0; i < n; ++i) 
	{
            xk[i] = X[IDX(i,k,ns)];
        }

        ukf->fx(xk, u, dt);
        for (int i = 0; i < n; ++i) 
	{
            ukf->sigma_x[IDX(i,k,ns)] = xk[i];
        }
    }

    /* 3. Compute predicted mean x */
    vec_zero(ukf->x, n);
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float w = ukf->Wm[k];
        for (int i = 0; i < n; ++i) 
	{
            ukf->x[i] += w * ukf->sigma_x[IDX(i,k,ns)];
        }
    }

    /* 4. Compute predicted covariance P */
    mat_zero(ukf->P, n, n);
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float wc = ukf->Wc[k];
        ukf_float dx[UKF_MAX_N];

        for (int i = 0; i < n; ++i) 
	{
            dx[i] = ukf->sigma_x[IDX(i,k,ns)] - ukf->x[i];
        }

        /* P += wc * dx * dx^T */
        for (int i = 0; i < n; ++i) 
	{
            for (int j = 0; j < n; ++j) 
	    {
                ukf->P[IDX(i,j,n)] += wc * dx[i] * dx[j];
            }
        }
    }

    /* Add process noise Q */
    for (int i = 0; i < n*n; ++i) 
    {
        ukf->P[i] += ukf->Q[i];
    }

    return UKF_OK;
}


/*!
 *----------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		ukf_status_t ukf_update(ukf_t *ukf, const ukf_float *z_meas)
 *
 *  @brief	Update step: incorporates measurement z using hx.
 *  
 *  @param	z: measurement vector of size n_z
 *
 *----------------------------------------------------------------------------------------------------------------------
 */
ukf_status_t ukf_update(ukf_t *ukf, const ukf_float *z_meas)
{
    if (!ukf || !ukf->hx || !z_meas) return UKF_ERR_DIM;

    int n = ukf->n_x;
    int m = ukf->n_z;
    int ns = ukf->n_sigma;

    /* 1. Transform sigma points into measurement space */
    ukf_float Z[UKF_MAX_M * UKF_MAX_SIGMA];
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float xk[UKF_MAX_N];
        for (int i = 0; i < n; ++i) 
	{
            xk[i] = ukf->sigma_x[IDX(i,k,ns)];
        }

        ukf_float zk[UKF_MAX_M];
        ukf->hx(xk, zk);

        for (int j = 0; j < m; ++j) 
	{
            Z[IDX(j,k,ns)] = zk[j];
        }
    }

    /* 2. Predicted measurement mean z_pred */
    ukf_float z_pred[UKF_MAX_M];
    vec_zero(z_pred, m);
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float w = ukf->Wm[k];
        for (int j = 0; j < m; ++j) 
	{
            z_pred[j] += w * Z[IDX(j,k,ns)];
        }
    }

    /* 3. Measurement covariance S */
    ukf_float S[UKF_MAX_M * UKF_MAX_M];
    mat_zero(S, m, m);
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float wc = ukf->Wc[k];
        ukf_float dz[UKF_MAX_M];
        for (int j = 0; j < m; ++j) 
	{
            dz[j] = Z[IDX(j,k,ns)] - z_pred[j];
	}

        for (int i = 0; i < m; ++i) 
	{
            for (int j = 0; j < m; ++j) 
	    {
                S[IDX(i,j,m)] += wc * dz[i] * dz[j];
            }
        }
    }

    /* Add measurement noise R */
    for (int i = 0; i < m*m; ++i) 
    {
        S[i] += ukf->R[i];
    }

    /* 4. Cross covariance P_xz */
    ukf_float P_xz[UKF_MAX_N * UKF_MAX_M];
    mat_zero(P_xz, n, m);
    for (int k = 0; k < ns; ++k) 
    {
        ukf_float wc = ukf->Wc[k];
        ukf_float dx[UKF_MAX_N];
        ukf_float dz[UKF_MAX_M];

        for (int i = 0; i < n; ++i) 
	{
            dx[i] = ukf->sigma_x[IDX(i,k,ns)] - ukf->x[i];
        }

        for (int j = 0; j < m; ++j) 
	{
            dz[j] = Z[IDX(j,k,ns)] - z_pred[j];
        }

        for (int i = 0; i < n; ++i) 
	{
            for (int j = 0; j < m; ++j) 
	    {
                P_xz[IDX(i,j,m)] += wc * dx[i] * dz[j];
            }
        }
    }

    /* 5. Kalman gain K = P_xz * S^-1 */
    ukf_float S_inv[UKF_MAX_M * UKF_MAX_M];
    if (mat_inv(S, S_inv, m) != 0) 
    {
        return UKF_ERR_INV;
    }

    ukf_float K[UKF_MAX_N * UKF_MAX_M];
    mat_mul(K, P_xz, S_inv, n, m, m);

    /* 6. Update state x and covariance P */
    ukf_float y[UKF_MAX_M]; /* innovation */
    for (int j = 0; j < m; ++j) 
    {
        y[j] = z_meas[j] - z_pred[j];
    }

    /* x = x + K*y */
    for (int i = 0; i < n; ++i) 
    {
        ukf_float sum = 0.0;
        for (int j = 0; j < m; ++j) 
	{
            sum += K[IDX(i,j,m)] * y[j];
        }
        ukf->x[i] += sum;
    }

    /* P = P - K*S*K^T */
    ukf_float KS[UKF_MAX_N * UKF_MAX_M];
    ukf_float KSKT[UKF_MAX_N * UKF_MAX_N];

    mat_mul(KS, K, S, n, m, m);
    ukf_float KT[UKF_MAX_M * UKF_MAX_N];
    mat_transpose(KT, K, n, m);
    mat_mul(KSKT, KS, KT, n, m, n);

    for (int i = 0; i < n*n; ++i) 
    {
        ukf->P[i] -= KSKT[i];
    }

    return UKF_OK;
}

