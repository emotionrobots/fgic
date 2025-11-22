# ECM Battery Model with 20-Point Tables and Online Parameter Update

This module implements a **single-cell equivalent circuit model (ECM)** in C with:

- 20-element lookup tables for $OCV$, $H_\text{chg}$, $H_\text{dsg}$, $R_0$, $R_1$, $C_1$ vs $SOC$ (at $20^\circ\mathrm{C}$),
- **Arrhenius temperature compensation** for $R_0$, $R_1$, $C_1$,
- A **lumped thermal model** for cell temperature $T$,
- Online updates of:
  - $R_0$ from step response $dV/dI$ at rest entry,
  - $R_1$ and $C_1$ from a least-squares (LSQ) fit of the RC branch voltage $V_{RC}(t)$,
  - $H_\text{chg}$ and $H_\text{dsg}$ from UKF-estimated hysteresis.

The C API is designed to work with the **generic UKF module**:

- Direct access to $SOC$,
- Interpolated evaluation of $OCV(SOC, T)$, $R_0(SOC, T)$, $R_1(SOC, T)$, $C_1(SOC, T)$,
- Interpolated $H_\text{chg}(SOC)$, $H_\text{dsg}(SOC)$.

---

## 1. ECM Structure

The ECM uses:

- A **series resistance** $R_0(T, SOC)$,
- A **parallel RC branch** $(R_1(T, SOC), C_1(T, SOC))$,
- Open-circuit voltage $OCV(SOC)$,
- **Hysteresis** $H_\text{chg}(SOC)$ and $H_\text{dsg}(SOC)$,
- A **lumped thermal model** for core temperature.

The terminal voltage is:

$$
V = OCV(SOC) + H_\text{dir}(SOC) + V_{RC} - I R_0(SOC, T)
$$

where $I$ is current ($I > 0$ = discharge, $I < 0$ = charge), and $H_\text{dir}(SOC)$ is chosen from:

- $H_\text{chg}(SOC)$ when charging,
- $H_\text{dsg}(SOC)$ when discharging.

### 1.1 20-Point SOC Tables

Each of $OCV$, $H_\text{chg}$, $H_\text{dsg}$, $R_0$, $R_1$, $C_1$ is stored as a 20-element table sampled at:

$$
SOC_i = \frac{i}{19}, \quad i = 0, 1, \dots, 19
$$

At runtime, the module uses **linear interpolation** in $SOC$ to evaluate the functions between grid points.

---

## 2. State Dynamics

The dynamic state vector is:

$$
\mathbf{x} =
\begin{bmatrix}
SOC \\
V_{RC} \\
T
\end{bmatrix}
$$

where:

- $SOC \in [0,1]$ is state of charge,
- $V_{RC}$ is the RC branch voltage,
- $T$ is the cell core temperature in $\ ^\circ\mathrm{C}$.

### 2.1 SOC Dynamics

Let $Q_\text{cell}$ be the capacity in Ah and $Q_C = 3600 Q_\text{cell}$ in Coulombs. With the convention $I > 0$ as discharge:

$$
\frac{d(SOC)}{dt} = -\frac{I}{Q_C}
$$

Discretized with step $\Delta t$:

$$
SOC_{k+1} = SOC_k - \frac{\Delta t}{Q_C} I_k
$$

The code clamps $SOC$ to $[0,1]$.

### 2.2 RC Branch Dynamics

The RC branch obeys:

$$
\frac{dV_{RC}}{dt} = -\frac{1}{R_1 C_1} V_{RC} + \frac{1}{C_1} I
$$

Using forward Euler discretization:

$$
V_{RC, k+1} = V_{RC,k}
+ \Delta t \left(
-\frac{V_{RC,k}}{R_1 C_1} + \frac{I_k}{C_1}
\right)
$$

Here $R_1$ and $C_1$ are computed from the tables and temperature (see Arrhenius scaling below).

### 2.3 Thermal Dynamics

The lumped thermal model uses:

$$
\frac{dT}{dt} = \frac{1}{C_{th}} \left( I^2 R_0(T, SOC) - \frac{T - T_{amb}}{R_{th}} \right)
$$

where:

- $C_{th}$ is the thermal capacitance (J/°C),
- $R_{th}$ is a thermal resistance (°C/W),
- $T_{amb}$ is ambient temperature.

Discretization:

$$
T_{k+1} = T_k + \Delta t \,
\frac{1}{C_{th}}
\left(
I_k^2 R_0(T_k, SOC_k) - \frac{T_k - T_{amb,k}}{R_{th}}
\right)
$$

---

## 3. Arrhenius Temperature Compensation

The tables for $R_0$, $R_1$, $C_1$ are defined at a reference temperature $T_\text{ref} = 20^\circ\mathrm{C}$. For a generic parameter $k(T)$ (e.g. $R_0$ at a given SOC), we use:

$$
k(T) = k_\text{ref} \exp\left(
-\frac{E_a}{R_g} \left(
\frac{1}{T_K} - \frac{1}{T_{\text{ref},K}}
\right)
\right)
$$

where:

- $k_\text{ref}$ is the value from the table at $T_\text{ref}$,
- $E_a$ is the activation energy (J/mol),
- $R_g$ is the gas constant,
- $T_K = T + 273.15$ and $T_{\text{ref},K} = T_\text{ref} + 273.15$.

This is applied separately for $R_0$, $R_1$, and $C_1$:

- $R_0(SOC, T)$: from $R_0$ table and $E_{a,R0}$,
- $R_1(SOC, T)$: from $R_1$ table and $E_{a,R1}$,
- $C_1(SOC, T)$: from $C_1$ table and $E_{a,C1}$.

In code, each table entry is interpolated in $SOC$ to obtain $k_\text{ref}$, then scaled by the Arrhenius factor.

---

## 4. Hysteresis Model

The model uses static *charge* and *discharge* hysteresis tables:

- $H_\text{chg}(SOC)$,
- $H_\text{dsg}(SOC)$.

The active hysteresis term $H_\text{dir}(SOC)$ is chosen according to current direction:

- If $I > 0$ (discharge), use $H_\text{dsg}(SOC)$,
- If $I < 0$ (charge), use $H_\text{chg}(SOC)$,
- If $I \approx 0$, the code keeps the last direction (charge or discharge).

The terminal voltage equation is:

$$
V = OCV(SOC) + H_\text{dir}(SOC) + V_{RC} - I R_0(SOC, T)
$$

### 4.1 Updating Hysteresis from UKF

If a UKF maintains a dynamic hysteresis state $H_\text{UKF}$, the ECM can **calibrate** its static tables from this estimate:

- Given a hysteresis estimate $H_\text{est}$ at some $SOC$ and a flag indicating charge vs discharge,
- The nearest SOC bin is updated using an exponential moving average:

$$
H_{\text{new}} = (1 - \alpha) H_{\text{old}} + \alpha H_\text{est}
$$

This is implemented by `ecm_update_h_from_ukf()` for both $H_\text{chg}$ and $H_\text{dsg}$.

---

## 5. Online Parameter Updates

The ECM supports **online identification** of $R_0$, $R_1$, and $C_1$ from operating data.

### 5.1 $R_0$ from $dV/dI$ at Rest Entry

When current steps from a load $I_\text{load}$ to rest ($I \approx 0$):

- The RC branch voltage $V_{RC}$ cannot change instantaneously,
- The terminal voltage exhibits an immediate step of approximately $I_\text{load} R_0$.

Let $V_{\text{prev}}$ be the measured voltage at the last load step, and $V_{\text{rest}}$ be the first measurement at rest. The code approximates:

$$
\Delta V = V_{\text{rest}} - V_{\text{prev}}, \quad
\Delta I = 0 - I_\text{load}
$$

Then:

$$
R_0 \approx -\frac{\Delta V}{\Delta I}
$$

This estimate is assigned to the nearest SOC bin using an exponential smoother:

$$
R_{0,\text{new}} = (1 - \alpha) R_{0,\text{old}} + \alpha R_{0,\text{est}}
$$

The logic is triggered in `ecm_update_from_measurement()` when the code detects a transition from $|I| > I_\text{thr}$ to $|I| \le I_\text{thr}$.

### 5.2 $R_1$ and $C_1$ from LSQ Fit of $V_{RC}(t)$

After the current step to rest ($I \approx 0$), the RC branch voltage decays exponentially:

$$
V_{RC}(t) = V_{RC}(0^+) e^{-t/\tau}, \quad \tau = R_1 C_1
$$

Taking the absolute value and logarithm:

$$
\ln |V_{RC}(t)| = \ln V_{RC}(0^+) - \frac{1}{\tau} t
$$

The code collects samples $\{t_i, |V_{RC}(t_i)|\}$ during the rest interval and fits:

$$
y_i = \ln |V_{RC}(t_i)| \approx a + b t_i
$$

via least-squares (LSQ) regression. The slope $b$ gives:

$$
\tau = -\frac{1}{b}
$$

At the moment of entering rest, we also record:

- The pre-rest current $I_\text{step}$,
- The initial RC voltage $V_{RC}(0^+)$.

From the first-order RC model:

$$
V_{RC}(0^+) \approx I_\text{step} R_1
$$

so we can estimate:

$$
R_1 \approx \frac{V_{RC}(0^+)}{|I_\text{step}|}, \qquad
C_1 \approx \frac{\tau}{R_1}
$$

Just like $R_0$, the code updates the nearest SOC bin with a smoothing factor $\alpha$:

$$
R_{1,\text{new}} = (1 - \alpha) R_{1,\text{old}} + \alpha R_{1,\text{est}}
$$

$$
C_{1,\text{new}} = (1 - \alpha) C_{1,\text{old}} + \alpha C_{1,\text{est}}
$$

This behavior is all contained in `ecm_update_from_measurement()`.

---

## 6. C API Summary

Key APIs in `ecm.h`:

```c
/* Initialization and state */
void   ecm_init_default(ecm_t *ecm);
void   ecm_reset_state(ecm_t *ecm, double soc0, double T0);

/* Dynamics and terminal voltage */
void   ecm_step(ecm_t *ecm, double I, double T_amb, double dt);
double ecm_terminal_voltage(const ecm_t *ecm, double I);

/* Lookups for UKF models */
double ecm_get_soc(const ecm_t *ecm);
double ecm_lookup_ocv(const ecm_t *ecm, double soc, double T);
double ecm_lookup_r0 (const ecm_t *ecm, double soc, double T);
double ecm_lookup_r1 (const ecm_t *ecm, double soc, double T);
double ecm_lookup_c1 (const ecm_t *ecm, double soc, double T);
double ecm_lookup_h_chg(const ecm_t *ecm, double soc);
double ecm_lookup_h_dsg(const ecm_t *ecm, double soc);

/* Online updates */
void ecm_update_from_measurement(ecm_t *ecm,
                                 double I,
                                 double V_meas,
                                 double vrc_est,
                                 double dt);

void ecm_update_h_from_ukf(ecm_t *ecm,
                           double soc,
                           double H_est,
                           int is_chg);
```

---

# Explanation of `ecm_test.c`

This document explains the structure and behavior of the **ECM unit test** implemented in `ecm_test.c`.  
The test exercises:

- 20-point table–based ECM,
- Arrhenius temperature compensation,
- Online parameter updates for $R_0$, $R_1$, and $C_1$,
- Behavior across **charge → rest** and **discharge → rest** transitions.

---

## 1. Purpose of `ecm_test.c`

The file `ecm_test.c` builds a **controlled experiment** to validate:

1. That the **ECM dynamics** (SOC, $V_{RC}$, $T$) behave sensibly under charge/discharge and rest.
2. That the **online parameter update mechanisms** can adjust:
   - $R_0$ from voltage step responses at rest,
   - $R_1$ and $C_1$ from the decay of the RC branch voltage $V_{RC}(t)$.
3. That a **model ECM** can be driven toward a **true ECM** using only:
   - Measured terminal voltage $V_\text{meas}$,
   - Estimated $V_{RC}$ (e.g., from a UKF; in the test, we “cheat” and use the true $V_{RC}$ to isolate the identification logic).

---

## 2. Two ECM Instances: True vs Model

`ecm_test.c` creates two ECM objects:

- `e_true` — the **true cell**:
  - Starts from `ecm_init_default()`,
  - Then its tables for $R_0$, $R_1$, and $C_1$ are deliberately modified:
    - E.g. $R_0$ is made larger, $R_1$ smaller, $C_1$ larger than the model.
  - This represents the “real” battery.

- `e_model` — the **model ECM**:
  - Uses the default table values,
  - Represents the **initial guess** used for modeling and identification.

Both are initialized with:

```c
ecm_init_default(&e_true);
ecm_init_default(&e_model);

/* Modify true parameters so there's something to learn */
for (int i = 0; i < ECM_TABLE_SIZE; ++i) {
    e_true.r0_table[i] *= 1.5;
    e_true.r1_table[i] *= 0.7;
    e_true.c1_table[i] *= 1.4;
}
```

Dynamic states are then reset with different initial SOCs, e.g.:

e_true starts at $SOC = 0.8$,

e_model starts at $SOC = 0.7$.

---

## 3. Test Phases: Charge/Rest/Discharge/Rest

The test runs in four phases, each with a fixed current:

Phase 0 – Charging
$I = -1.0\ \mathrm{A}$ (charging, SOC increases).

Phase 1 – Rest after Charge
$I = 0.0\ \mathrm{A}$ (no current, $V_{RC}$ decays, useful for $R_1$, $C_1$ identification).

Phase 2 – Discharging
$I = +1.0\ \mathrm{A}$ (discharge, SOC decreases).

Phase 3 – Rest after Discharge
$I = 0.0\ \mathrm{A}$ (again a rest interval to identify parameters).

Each phase runs for a fixed number of steps (e.g. 100 samples per phase), with constant step size:

   $$ \delta t = 1 s $$

At each time step $k$:

Time is set to:

   $$ t_k = (k+1)\delta t $$

The true ECM and model ECM are both stepped forward using:

```c
ecm_step(&e_true,  I, T_amb, dt);
ecm_step(&e_model, I, T_amb, dt);
```

where $T_\text{amb}$ is the ambient temperature (fixed, e.g. $20^\circ\mathrm{C}$).

---

## 4. True vs Measured vs Model Voltage

After stepping the true ECM, the test computes the clean true terminal voltage:
```c
double V_true_clean = ecm_terminal_voltage(&e_true, I);
```

To simulate real measurements, the test then adds Gaussian noise:

   $$ V_text{meas} = V_\text{true} + n_V $$

where $n_V$ is normally distributed with, for example, $\sigma \approx 5\ \mathrm{mV}$.

In code:

```c
double V_meas = V_true_clean + rand_normal(0.005);
```

The model ECM then computes its own terminal voltage:

```c
double V_model = ecm_terminal_voltage(&e_model, I);
```

---

## . Online Parameter Updating

The crux of ecm_test.c is calling:

```c
ecm_update_from_measurement(&e_model, I, V_meas, vrc_est, dt);
```

where:

* I is the current at this step,

* V_meas is the noisy terminal voltage from the true ECM,

* vrc_est is an estimate of the RC branch voltage $V_{RC}$,

* dt is the time step.

In the test, to isolate the parameter identification logic, vrc_est is taken from the true model’s internal state:

```c
double vrc_est = e_true.v_rc;  /* As if from a perfect UKF */
```

This function performs three main tasks:

### 5.1 Detecting Rest Entry and Updating $R_0$

When the current transitions from a non-rest region ($|I| > I_\text{rest}$) to a rest region ($|I| \le I_\text{rest}$), the terminal voltage shows an instantaneous step approximately equal to:

   $$ \triangle_V \approx -I_\text{load} R_0 $$


