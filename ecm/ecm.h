/*!
 *=====================================================================================================================
 *
 * @file		ecm.h
 *
 * @brief		ECM header file
 *
 *=====================================================================================================================
 */
#ifndef __ECM_H__
#define __ECM_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ECM_TABLE_SIZE 		20
#define ECM_LSQ_MAX    		64   /* max samples for LSQ on VRC decay */

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Lumped ECM model with lookup tables and Arrhenius temperature scaling. 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
typedef struct {
    /* SOC grid (monotonic, usually 0..1) */
    double soc_grid[ECM_TABLE_SIZE];

    /* Base tables defined at reference temperature T_ref_C (20 °C) */
    double ocv_table[ECM_TABLE_SIZE];      /* OCV(SOC, T_ref) */
    double h_chg_table[ECM_TABLE_SIZE];    /* H_chg(SOC) at T_ref */
    double h_dsg_table[ECM_TABLE_SIZE];    /* H_dsg(SOC) at T_ref */
    double r0_table[ECM_TABLE_SIZE];       /* R0(SOC, T_ref) */
    double r1_table[ECM_TABLE_SIZE];       /* R1(SOC, T_ref) */
    double c1_table[ECM_TABLE_SIZE];       /* C1(SOC, T_ref) */

    /* Arrhenius activation energies (J/mol) for R0, R1, C1 */
    double Ea_R0;
    double Ea_R1;
    double Ea_C1;

    /* Reference temperature for tables (°C), typically 20 °C */
    double T_ref_C;

    /* Cell capacity (Ah) */
    double Q_Ah;

    /* Thermal parameters */
    double C_th;    /* thermal capacitance [J/°C] */
    double R_th;    /* thermal resistance [°C/W] */

    /* Dynamic states */
    double soc;     /* current SOC (0..1) */
    double v_rc;    /* RC branch voltage [V] */
    double T;       /* core cell temperature [°C] */

    /* Hysteresis direction memory: -1 charge, +1 discharge, 0 rest */
    int last_dir;

    /* Quit current that determine REST */
    double quit_current;

    /* --- Online parameter update bookkeeping --- */

    /* For R0 update via dV/dI at rest entry */
    double prev_I;
    double prev_V;
    int    prev_is_rest;

    /* For R1/C1 LSQ on VRC decay */
    int    in_rest_segment;
    double rest_time;                          /* time since rest entry [s] */
    double step_I;                             /* current before entering rest */
    double vrc0;                               /* |VRC| at start of rest */
    double t_hist[ECM_LSQ_MAX];                /* times within rest */
    double vrc_hist[ECM_LSQ_MAX];              /* |VRC| samples */
    int    hist_len;
} ecm_t;


/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Initialize ECM with default example tables (20 points, 0..1 SOC).
 * Sets reasonable defaults for OCV, R0, R1, C1, hysteresis, Arrhenius,
 * thermal parameters, and initial state.
 * Initialize ECM with default example tables (20 points, 0..1 SOC).
 *--------------------------------------------------------------------------------------------------------------------- 
 */
void ecm_init_default(ecm_t *ecm);


/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Reset dynamic state only (SOC, VRC, T, direction and update buffers). 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
void ecm_reset_state(ecm_t *ecm, double soc0, double T0);


/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Advance ECM dynamics for one time step.
 *
 * Inputs:
 *  ecm   : ECM instance
 *  I     : cell current [A], I > 0 = discharge, I < 0 = charge
 *  T_amb : ambient temperature [°C]
 *  dt    : time step [s]
 *
 * States updated: soc, v_rc, T, last_dir
 * Parameters R0/R1/C1 are read via lookup + Arrhenius scaling.
 *--------------------------------------------------------------------------------------------------------------------- 
 */
void ecm_step(ecm_t *ecm, double I, double T_amb, double dt);

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Compute terminal voltage for given current I using current ECM state. 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_terminal_voltage(const ecm_t *ecm, double I);



/* --- Lookup helpers for use in UKF models --- */

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Return internal SOC state. 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_get_soc(const ecm_t *ecm);


/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Interpolated lookups at arbitrary SOC, T (°C).
 * OCV is assumed weakly temperature dependent, so T is ignored here,
 * but passed for a consistent API.
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_lookup_ocv(const ecm_t *ecm, double soc, double T);


/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * R0, R1, C1 with Arrhenius temperature scaling from T_ref_C to T (°C). 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_lookup_r0(const ecm_t *ecm, double soc, double T);
double ecm_lookup_r1(const ecm_t *ecm, double soc, double T);
double ecm_lookup_c1(const ecm_t *ecm, double soc, double T);

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Hysteresis tables (no temperature scaling here). 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_lookup_h_chg(const ecm_t *ecm, double soc);
double ecm_lookup_h_dsg(const ecm_t *ecm, double soc);

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Convenience getters for current internal state: 
 *--------------------------------------------------------------------------------------------------------------------- 
 */
double ecm_get_ocv_now(const ecm_t *ecm);
double ecm_get_r0_now(const ecm_t *ecm);
double ecm_get_r1_now(const ecm_t *ecm);
double ecm_get_c1_now(const ecm_t *ecm);

/* --- Online parameter update APIs -------------------------------------- */

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Update internal R0, R1, C1 tables from measurements.
 *
 * Inputs:
 *  ecm      : ECM instance (the "model" you want to adapt)
 *  I        : current [A]
 *  V_meas   : measured terminal voltage [V]
 *  vrc_est  : estimated or measured RC-branch voltage [V] (e.g. from UKF),
 *             used for fitting R1/C1 during rest segments.
 *  dt       : time step [s]
 *
 * Effects:
 *  - On entering rest (|I| small) after a load current: update R0 SOC bin
 *    using dV/dI.
 *  - During rest, collect VRC(t) samples.
 *  - On exiting rest: LSQ fit of exponential decay VRC(t) to update
 *    R1 and C1 SOC bins via tau = R1*C1 and VRC(0+) = I_step * R1.
 *--------------------------------------------------------------------------------------------------------------------- 
 */
void ecm_update_from_measurement(ecm_t *ecm,
                                 double I,
                                 double V_meas,
                                 double vrc_est,
                                 double dt);

/*!
 *--------------------------------------------------------------------------------------------------------------------- 
 * Update H_chg / H_dsg tables from UKF-estimated hysteresis.
 *
 * Inputs:
 *  soc       : SOC at which hysteresis estimate applies
 *  H_est     : hysteresis voltage estimate [V]
 *  is_chg    : non-zero => charging hysteresis (H_chg), else discharging (H_dsg)
 *
 * A simple exponential moving average is applied at the nearest SOC bin.
 *--------------------------------------------------------------------------------------------------------------------- 
 */
void ecm_update_h_from_ukf(ecm_t *ecm,
                           double soc,
                           double H_est,
                           int is_chg);

#ifdef __cplusplus
}
#endif

#endif /* __ECM_H__ */

