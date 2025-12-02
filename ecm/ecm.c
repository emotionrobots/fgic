/*!
 *=====================================================================================================================
 *
 * @file		ecm.c
 *
 * @brief		Equivalent circuit model
 *
 *=====================================================================================================================
 */
#include "ecm.h"
#include <math.h>
#include <string.h>


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double clamp(double x, double lo, double hi)
 *
 * @brief	Clamping function 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
double clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double interp_soc(const double *grid, const double *tbl, int n, double soc)
 *
 * @brief	Linear interpolation on SOC grid (assumed monotonic). 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
double interp_soc(const double *grid, const double *tbl, int n, double soc)
{
    soc = clamp(soc, grid[0], grid[n - 1]);

    /* Below first bin? */
    if (soc <= grid[0]) return tbl[0];
    if (soc >= grid[n - 1]) return tbl[n - 1];

    for (int i = 0; i < n - 1; ++i) {
        double s0 = grid[i];
        double s1 = grid[i + 1];
        if (soc >= s0 && soc <= s1) {
            double t = (soc - s0) / (s1 - s0);
            return tbl[i] + t * (tbl[i + 1] - tbl[i]);
        }
    }
    /* Fallback (should not happen): */
    return tbl[n - 1];
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		int nearest_soc_idx(const double *grid, int n, double soc)
 *
 * @brief	Find nearest SOC bin index. 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
int nearest_soc_idx(const double *grid, int n, double soc)
{
    soc = clamp(soc, grid[0], grid[n - 1]);
    int best = 0;
    double best_err = fabs(soc - grid[0]);
    for (int i = 1; i < n; ++i) {
        double e = fabs(soc - grid[i]);
        if (e < best_err) {
            best_err = e;
            best = i;
        }
    }
    return best;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double arrhenius_scale(double k_ref, double Ea, double T_C, double Tref_C)
 *
 * @brief	Arrhenius scaling k(T) = k_ref * exp(-Ea/R * (1/T - 1/T_ref)).  T, T_ref in °C; internal 
 *              conversion to Kelvin.
 *
 * @param	k_ref 		Reference value at T=T_ref_C
 * @param	Ea		Energy coefficient; Ea > 0 k will grow with temperature; 
 *                                                  Ea < 0 k will shrink with temperature
 *
 * Rg		Thermal resistance TODO: check this value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
double arrhenius_scale(double k_ref, double Ea, double T_C, double Tref_C)
{
    /* If Ea ~ 0, skip scaling. */
    // if (fabs(Ea) < 1.0) return k_ref;

    const double Rg = 8.314462618; /* J/mol/K */
    double T  = T_C     + 273.15;
    double Tr = Tref_C  + 273.15;

    if (T < 1.0)  T  = 1.0;		/* Dont accept near 0 deg K */
    if (Tr < 1.0) Tr = 1.0;            

    double factor = exp( -Ea / Rg * (1.0 / T - 1.0 / Tr) );
    return k_ref * factor;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn  	void ecm_init_default(ecm_t *ecm)
 *
 * @brief	Init ECM model to default values
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
void ecm_init_default(ecm_t *ecm)
{
    memset(ecm, 0, sizeof(*ecm));

    /* Default SOC grid: 0..1 */
    for (int i = 0; i < ECM_TABLE_SIZE; ++i) {
        ecm->soc_grid[i] = (double)i / (double)(ECM_TABLE_SIZE - 1);
    }

    /* Example OCV / hysteresis / R0 / R1 / C1 tables at 20°C. */
    /* TODO: replace with actual tbl values */
    for (int i = 0; i < ECM_TABLE_SIZE; ++i) 
    {
        double s = ecm->soc_grid[i];

        /* Simple linear-ish OCV: 3.0 V at 0% -> 3.7 V at 100% */
        ecm->ocv_table[i] = 3.0 + 0.7 * s;

        /* Hysteresis magnitude a bit larger at low SOC (pure example) */
        ecm->h_chg_table[i] =  0.03 * (1.0 - s); /* charge branch */
        ecm->h_dsg_table[i] = -0.03 * (1.0 - s); /* discharge branch */

        /* R0 decreases slightly with SOC */
        ecm->r0_table[i] = 0.030 - 0.010 * s;        /* ohm */

        /* R1 decreases with SOC, C1 increases a bit */
        ecm->r1_table[i] = 0.015 - 0.008 * s;       /* ohm */
        ecm->c1_table[i] = 800.0 + 400.0 * s;       /* F */
    }

    /* Arrhenius parameters (example values) */
    ecm->Ea_R0   = -20000.0;  // shrinking with temp 
    ecm->Ea_R1   = -15000.0;  // shrinking with temp
    ecm->Ea_C1   =  10000.0;  // growing with temp

    /* Reference temperature for tables */
    ecm->T_ref_C = 25.0;

    /* Capacity and thermal */
    ecm->Q_Ah = 2.5;      /* 2.5 Ah cell (example) */
    ecm->C_th = 200.0;    /* J/°C */
    ecm->R_th = 3.0;      /* °C/W */

    /* Quit current */
    ecm->quit_current = 1e-3;   /* 1mA */

    /* Dynamic state defaults */
    ecm->soc = 0.5;
    ecm->v_rc = 0.0;
    ecm->T   = 20.0;
    ecm->last_dir = 0;

    ecm->prev_I = 0.0;
    ecm->prev_V = ecm->ocv_table[0];
    ecm->prev_is_rest = 1;

    ecm->in_rest_segment = 0;
    ecm->rest_time = 0.0;
    ecm->step_I = 0.0;
    ecm->vrc0 = 0.0;
    ecm->hist_len = 0;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		void ecm_reset_state(ecm_t *ecm, double soc0, double T0)
 *
 * @brief	Reset ECM state but not model coef.
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
void ecm_reset_state(ecm_t *ecm, double soc0, double T0)
{
    ecm->soc = clamp(soc0, 0.0, 1.0);
    ecm->v_rc = 0.0;
    ecm->T = T0;
    ecm->last_dir = 0;

    ecm->prev_I = 0.0;
    ecm->prev_V = ecm_lookup_ocv(ecm, ecm->soc, ecm->T);
    ecm->prev_is_rest = 1;

    ecm->in_rest_segment = 0;
    ecm->rest_time = 0.0;
    ecm->step_I = 0.0;
    ecm->vrc0 = 0.0;
    ecm->hist_len = 0;
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_get_soc(const ecm_t *ecm)
 *
 * @brief	Get ECM SOC
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_get_soc(const ecm_t *ecm)
{
    return ecm->soc;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_ocv(const ecm_t *ecm, double soc, double T)
 *
 * @brief	Lookup OCV given SOC and T
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_ocv(const ecm_t *ecm, double soc, double T)
{
    (void)T; /* ignoring temperature for OCV in this example */
    return interp_soc(ecm->soc_grid, ecm->ocv_table, ECM_TABLE_SIZE, soc);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_h_chg(const ecm_t *ecm, double soc)
 * 
 * @brief	Lookup H_chg table value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_h_chg(const ecm_t *ecm, double soc)
{
    return interp_soc(ecm->soc_grid, ecm->h_chg_table, ECM_TABLE_SIZE, soc);
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_h_dsg(const ecm_t *ecm, double soc)
 * 
 * @brief	Lookup H_dsg table value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_h_dsg(const ecm_t *ecm, double soc)
{
    return interp_soc(ecm->soc_grid, ecm->h_dsg_table, ECM_TABLE_SIZE, soc);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_r0(const ecm_t *ecm, double soc)
 * 
 * @brief	Lookup R0 table value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_r0(const ecm_t *ecm, double soc, double T)
{
    double r0_ref = interp_soc(ecm->soc_grid, ecm->r0_table, ECM_TABLE_SIZE, soc);
    return arrhenius_scale(r0_ref, ecm->Ea_R0, T, ecm->T_ref_C);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_r1(const ecm_t *ecm, double soc)
 * 
 * @brief	Lookup R1 table value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_r1(const ecm_t *ecm, double soc, double T)
{
    double r1_ref = interp_soc(ecm->soc_grid, ecm->r1_table, ECM_TABLE_SIZE, soc);
    return arrhenius_scale(r1_ref, ecm->Ea_R1, T, ecm->T_ref_C);
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		double ecm_lookup_c1(const ecm_t *ecm, double soc)
 * 
 * @brief	Lookup C1 table value
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_lookup_c1(const ecm_t *ecm, double soc, double T)
{
    double c1_ref = interp_soc(ecm->soc_grid, ecm->c1_table, ECM_TABLE_SIZE, soc);
    return arrhenius_scale(c1_ref, ecm->Ea_C1, T, ecm->T_ref_C);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double ecm_get_ocv_now(const ecm_t *ecm)
 *  
 *  @brief	Convenience getters at current state 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_get_ocv_now(const ecm_t *ecm)
{
    return ecm_lookup_ocv(ecm, ecm->soc, ecm->T);
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double ecm_get_r0_now(const ecm_t *ecm)
 *  
 *  @brief	Convenience getters at current state 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_get_r0_now(const ecm_t *ecm)
{
    return ecm_lookup_r0(ecm, ecm->soc, ecm->T);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double ecm_get_r1_now(const ecm_t *ecm)
 *  
 *  @brief	Convenience getters at current state 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_get_r1_now(const ecm_t *ecm)
{
    return ecm_lookup_r1(ecm, ecm->soc, ecm->T);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double ecm_get_c1_now(const ecm_t *ecm)
 *  
 *  @brief	Convenience getters at current state 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_get_c1_now(const ecm_t *ecm)
{
    return ecm_lookup_c1(ecm, ecm->soc, ecm->T);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void ecm_step(ecm_t *ecm, double I, double T_amb, double dt)
 *
 *  @brief	Compute the ECM dynamics; normally this require UKF
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
void ecm_step(ecm_t *ecm, double I, double T_amb, double dt)
{
    /* Capacity in Coulombs */
    double Qc = ecm->Q_Ah * 3600.0;

    /* SOC update: I>0 discharge, I<0 charge */
    double dSOC = -(I * dt) / Qc;
    ecm->soc = clamp(ecm->soc + dSOC, 0.0, 1.0);

    /* RC branch update */
    double R1 = ecm_lookup_r1(ecm, ecm->soc, ecm->T);
    double C1 = ecm_lookup_c1(ecm, ecm->soc, ecm->T);
    double tau1 = R1 * C1;

    if (tau1 < 1e-9) tau1 = 1e-9; /* avoid degenerate */

    double dVRC = dt * ( -ecm->v_rc / tau1 + I / C1 );
    ecm->v_rc += dVRC;

    /* Thermal update */
    double R0 = ecm_lookup_r0(ecm, ecm->soc, ecm->T);
    double power_loss = I * I * R0;
    double dT = dt * ( (power_loss - (ecm->T - T_amb) / ecm->R_th) / ecm->C_th );
    ecm->T += dT;

    /* Track direction for hysteresis sign */
    if (I > ecm->quit_current) 
    {
        ecm->last_dir = +1; /* discharge */
    } 
    else if (I < -ecm->quit_current) 
    {
        ecm->last_dir = -1; /* charge */
    }
    else
    {
        ecm->last_dir = 0; /* rest */
    }

}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double ecm_terminal_voltage(const ecm_t *ecm, double I)
 *
 *  @brief	Compute ECM terminal voltage
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
double ecm_terminal_voltage(const ecm_t *ecm, double I)
{
    double ocv = ecm_lookup_ocv(ecm, ecm->soc, ecm->T);

    /* Choose hysteresis branch based on direction or last_dir */
    const double I_eps = 1e-3;
    int dir = ecm->last_dir;
    if (I > I_eps)      dir = +1;
    else if (I < -I_eps) dir = -1;

    double H = 0.0;
    if (dir > 0) {
        H = ecm_lookup_h_dsg(ecm, ecm->soc);
    } else if (dir < 0) {
        H = ecm_lookup_h_chg(ecm, ecm->soc);
    }

    double R0 = ecm_lookup_r0(ecm, ecm->soc, ecm->T);

    double V = ocv + H + ecm->v_rc - I * R0;
    return V;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 * @fn		void ecm_finalize_lsq_segment(ecm_t *ecm)
 *
 * @brief	Perform online parameter updates 
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
void ecm_finalize_lsq_segment(ecm_t *ecm)
{
    if (!ecm->in_rest_segment) return;

    if (ecm->hist_len < 3) 
    {
        ecm->in_rest_segment = 0;
        ecm->hist_len = 0;
        return;
    }

    /* Fit ln(VRC) = a + b t => tau = -1/b */
    double n = (double)ecm->hist_len;
    double sum_t = 0.0, sum_y = 0.0, sum_tt = 0.0, sum_ty = 0.0;

    for (int i = 0; i < ecm->hist_len; ++i) 
    {
        double t = ecm->t_hist[i];
        double v = ecm->vrc_hist[i];
        if (v <= 0.0) continue;
        double y = log(v);

        sum_t  += t;
        sum_y  += y;
        sum_tt += t * t;
        sum_ty += t * y;
    }

    double denom = n * sum_tt - sum_t * sum_t;
    if (fabs(denom) < 1e-12) 
    {
        ecm->in_rest_segment = 0;
        ecm->hist_len = 0;
        return;
    }

    double b = (n * sum_ty - sum_t * sum_y) / denom;

    /* need negative slope */
    if (b >= -1e-6) 
    { 
        ecm->in_rest_segment = 0;
        ecm->hist_len = 0;
        return;
    }

    double tau = -1.0 / b;

    /* Need initial VRC and step current for R1/C1 identification */
    if (fabs(ecm->step_I) < 1e-3 || ecm->vrc0 <= 0.0) 
    {
        ecm->in_rest_segment = 0;
        ecm->hist_len = 0;
        return;
    }

    /* R1 ~ VRC(0+)/I_step, C1 ~ tau / R1 */
    double R1_est = ecm->vrc0 / fabs(ecm->step_I);
    double C1_est = tau / R1_est;

    if (R1_est <= 0.0 || C1_est <= 0.0) 
    {
        ecm->in_rest_segment = 0;
        ecm->hist_len = 0;
        return;
    }

    /* Update nearest SOC bin */
    int idx = nearest_soc_idx(ecm->soc_grid, ECM_TABLE_SIZE, ecm->soc);
    double alpha = 0.3; /* adaptation rate */

    ecm->r1_table[idx] = (1.0 - alpha) * ecm->r1_table[idx] + alpha * R1_est;
    ecm->c1_table[idx] = (1.0 - alpha) * ecm->c1_table[idx] + alpha * C1_est;

    ecm->in_rest_segment = 0;
    ecm->hist_len = 0;
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void ecm_update_from_measurement(ecm_t *ecm, double I, double V_meas, double vrc_est, double dt)
 *
 *  @brief	Perform ECM update from measurement
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
void ecm_update_from_measurement(ecm_t *ecm, double I, double V_meas, double vrc_est, double dt)
{
    const double I_rest_thr = 0.02; /* 20 mA threshold for "rest" */
    const double VRC_min    = 1e-5;

    int is_rest = (fabs(I) <= I_rest_thr);

    /* ----- R0 update at rest entry (dV/dI) ----- */
    if (is_rest && !ecm->prev_is_rest && fabs(ecm->prev_I) > I_rest_thr) 
    {
        double dI = I - ecm->prev_I;  /* should be -prev_I */
        double dV = V_meas - ecm->prev_V;

        if (fabs(dI) > 1e-6) 
	{
            double R0_est = -dV / dI; /* sign so that positive R0 */
            if (R0_est > 1e-4 && R0_est < 1.0) 
	    {
                int idx = nearest_soc_idx(ecm->soc_grid, ECM_TABLE_SIZE, ecm->soc);
                double alpha = 0.3;
                ecm->r0_table[idx] = (1.0 - alpha) * ecm->r0_table[idx] + alpha * R0_est;
            }
        }

        /* Start new LSQ segment for VRC decay */
        ecm->in_rest_segment = 1;
        ecm->rest_time = 0.0;
        ecm->hist_len = 0;
        ecm->step_I = ecm->prev_I;
        ecm->vrc0 = fabs(vrc_est);

        if (ecm->vrc0 > VRC_min && ecm->hist_len < ECM_LSQ_MAX) 
	{
            ecm->t_hist[0] = 0.0;
            ecm->vrc_hist[0] = ecm->vrc0;
            ecm->hist_len = 1;
        }
    }

    /* ----- LSQ accumulation within rest segment ----- */
    if (is_rest && ecm->in_rest_segment) 
    {
        ecm->rest_time += dt;
        double vrc_abs = fabs(vrc_est);

        if (vrc_abs > VRC_min && ecm->hist_len < ECM_LSQ_MAX) 
	{
            ecm->t_hist[ecm->hist_len] = ecm->rest_time;
            ecm->vrc_hist[ecm->hist_len] = vrc_abs;
            ecm->hist_len++;
        }
    }

    /* ----- Exiting rest: finalize LSQ for R1/C1 ----- */
    if (!is_rest && ecm->prev_is_rest && ecm->in_rest_segment) 
    {
        ecm_finalize_lsq_segment(ecm);
    }

    /* Save previous sample information */
    ecm->prev_I = I;
    ecm->prev_V = V_meas;
    ecm->prev_is_rest = is_rest;
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void ecm_update_h_from_ukf(ecm_t *ecm, double soc, double H_est, int is_chg)
 *  
 *  @brief	Update H from UKF
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
void ecm_update_h_from_ukf(ecm_t *ecm, double soc, double H_est, int is_chg)
{
    int idx = nearest_soc_idx(ecm->soc_grid, ECM_TABLE_SIZE, soc);
    double alpha = 0.2; /* hysteresis adaptation rate */

    if (is_chg) 
    {
        ecm->h_chg_table[idx] = (1.0 - alpha) * ecm->h_chg_table[idx] + alpha * H_est;
    } 
    else 
    {
        ecm->h_dsg_table[idx] = (1.0 - alpha) * ecm->h_dsg_table[idx] + alpha * H_est;
    }
}

