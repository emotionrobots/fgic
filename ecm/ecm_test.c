/*!
 *=====================================================================================================================
 *
 * @file		ecm_test.c
 *
 * @brief		ECM unit test
 *
 *=====================================================================================================================
 */
#include <stdio.h>
#include <math.h>
#include "ecm.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* Simple normal noise generator for measurements */
static unsigned int lcg_state = 123456789u;


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double rand_uniform(void)
 *
 *  @brief	Create random number with uniform distribution
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
double rand_uniform(void)
{
    lcg_state = 1664525u * lcg_state + 1013904223u;
    return (double)(lcg_state & 0xFFFFFFu) / (double)0xFFFFFFu;
}


/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		double rand_normal(double std)
 *
 *  @brief 	Create random number with normal distribution
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
double rand_normal(double std)
{
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    if (u1 < 1e-12) u1 = 1e-12;
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    return std * r * cos(theta);
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void print_header(void)
 *
 *  @brief	Printer header of CSV file
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
void print_header(void)
{
    printf("phase,k,t,I,"
           "soc_true,soc_model,"
           "V_true,V_model,"
           "R0_true,R0_model,"
           "R1_true,R1_model,"
           "C1_true,C1_model\n");
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *
 *  @fn		void run_phase(int phase_id, double I, int steps,
 *             		       ecm_t *e_true, ecm_t *e_model,
 *             		       double T_amb, double dt, int *pk)
 *
 *  @brief	Run simulation for a number of steps (phase) at given I, T_amb
 *
 *---------------------------------------------------------------------------------------------------------------------
 */
static 
void run_phase(int phase_id, double I, int steps,
               ecm_t *e_true, ecm_t *e_model,
               double T_amb, double dt, int *pk)
{
    for (int n = 0; n < steps; ++n, ++(*pk)) 
    {
        int k = *pk;
        double t = (k + 1) * dt;

	// Compute true ECM dynamic one step and find terminal voltage with 5mV noise added 
        ecm_step(e_true, I, T_amb, dt);
        double V_true_clean = ecm_terminal_voltage(e_true, I);
        double V_meas = V_true_clean + rand_normal(0.005);

	// Compute model ECM dynamic one step, update ECM model, and find model terminal voltage 
        ecm_step(e_model, I, T_amb, dt);
        double vrc_est = e_true->v_rc;
        ecm_update_from_measurement(e_model, I, V_meas, vrc_est, dt);
        double V_model = ecm_terminal_voltage(e_model, I);


        double R0_true  = ecm_lookup_r0(e_true,  e_true->soc,  e_true->T);
        double R0_model = ecm_lookup_r0(e_model, e_model->soc, e_model->T);
        double R1_true  = ecm_lookup_r1(e_true,  e_true->soc,  e_true->T);
        double R1_model = ecm_lookup_r1(e_model, e_model->soc, e_model->T);
        double C1_true  = ecm_lookup_c1(e_true,  e_true->soc,  e_true->T);
        double C1_model = ecm_lookup_c1(e_model, e_model->soc, e_model->T);

        printf("%d,%d,%.1f,%.4f,"
               "%.4f,%.4f,"
               "%.4f,%.4f,"
               "%.6f,%.6f,"
               "%.6f,%.6f,"
               "%.2f,%.2f\n",
               phase_id, k, t, I,
               e_true->soc, e_model->soc,
               V_true_clean, V_model,
               R0_true, R0_model,
               R1_true, R1_model,
               C1_true, C1_model);
    }
}

/*!
 *---------------------------------------------------------------------------------------------------------------------
 *  Main
 *---------------------------------------------------------------------------------------------------------------------
 */
int main(void)
{
    ecm_t e_true, e_model;
    ecm_init_default(&e_true);
    ecm_init_default(&e_model);

    /* ... same modifications as above ... */

    double dt = 1.0;
    double T_amb = 20.0;
    int k = 0;

    print_header();

    // Charge @1A for 100 steps
    run_phase(0, -1.0, 100, &e_true, &e_model, T_amb, dt, &k);

    // Rest for 100 steps 
    run_phase(1,  0.0, 100, &e_true, &e_model, T_amb, dt, &k);

    // Discharge @1A for 100 steps 
    run_phase(2,  1.0, 100, &e_true, &e_model, T_amb, dt, &k);

    // Rest 100 steps 
    run_phase(3,  0.0, 100, &e_true, &e_model, T_amb, dt, &k);

    return 0;
}

