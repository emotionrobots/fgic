/*!
 *======================================================================================================================
 *
 * @file	ecm.h
 *
 * @brief	ECM module header
 *
 *======================================================================================================================
 */
#ifndef __ECM_H__
#define __ECM_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ECM_TABLE_SIZE 			20
#define ECM_LSQ_MAX    			64   	/* max samples for LSQ on VRC decay */
#define ECM_QUIT_CURRENT		0.001   /* 1 mA */


/* Lumped ECM model with lookup tables and Arrhenius temperature scaling. */
typedef struct 
{
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

    /* Hysteresis direction memory: -1 charge, +1 discharge, 0 unknown */
    int last_dir;

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
    double I_quit;			       /* Quit current */
} 
ecm_t;


void ecm_init_default(ecm_t *ecm, double T0);
void ecm_reset_state(ecm_t *ecm, double soc0, double T0);
void ecm_step(ecm_t *ecm, double I, double T_amb, double dt);
double ecm_terminal_voltage(const ecm_t *ecm, double I);
double ecm_get_soc(const ecm_t *ecm);

double ecm_lookup_ocv(const ecm_t *ecm, double soc, double T);
double ecm_lookup_r0(const ecm_t *ecm, double soc, double T);
double ecm_lookup_r1(const ecm_t *ecm, double soc, double T);
double ecm_lookup_c1(const ecm_t *ecm, double soc, double T);

double ecm_lookup_h_chg(const ecm_t *ecm, double soc);
double ecm_lookup_h_dsg(const ecm_t *ecm, double soc);

double ecm_get_ocv_now(const ecm_t *ecm);
double ecm_get_r0_now(const ecm_t *ecm);
double ecm_get_r1_now(const ecm_t *ecm);
double ecm_get_c1_now(const ecm_t *ecm);

void ecm_update_from_measurement(ecm_t *ecm, double I, double V_meas, double vrc_est, double dt);
void ecm_update_h_from_ukf(ecm_t *ecm, double soc, double H_est, int is_chg);

#ifdef __cplusplus
}
#endif

#endif /* __ECM_H__ */

