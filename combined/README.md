# Combined UKF + ECM Battery Model Test

This document explains the **combined Unscented Kalman Filter (UKF) + Equivalent Circuit Model (ECM)** unit test implemented in `ukf_ecm_test.c`.

The goal is to show that:

- The **UKF** tracks the internal battery states:
  - State of charge $SOC$,
  - RC branch voltage $V_{RC}$,
  - Hysteresis $H$,
  - Core temperature $T$,
- The **ECM** parameter tables and Arrhenius scaling are used consistently,
- The **UKF+ECM measurement model** reproduces the terminal voltage of a “true” cell under a discharging–rest–charging–rest scenario.

---

## 1. Model Overview

### 1.1 UKF State and Measurement

The UKF state is a 4-dimensional vector:

$$ \mathbf{x} = 
\begin{bmatrix}
SOC \\
V_{RC} \\
H \\
T
\end{bmatrix}
$$

where:

- $SOC$ is the state of charge,
- $V_{RC}$ is the RC branch voltage across the diffusion branch,
- $H$ is the hysteresis voltage,
- $T$ is the core cell temperature in $\ ^\circ\mathrm{C}$.

The measurement is the **terminal voltage**:

$$
\mathbf{z} = [ V ].
$$

The UKF input is a 2-vector:

$$
\mathbf{u} =
\begin{bmatrix}
I \\
T_{amb}
\end{bmatrix},
$$

where $I$ is the current ($I>0$ = discharge, $I<0$ = charge) and $T_{amb}$ is the ambient temperature.

---

### 1.2 ECM Parameterization

The ECM uses 20-point tables vs $SOC$ for:

- $OCV(SOC)$,
- $H_{\text{chg}}(SOC)$ and $H_{\text{dsg}}(SOC)$,
- $R_0(SOC, T_\text{ref})$,
- $R_1(SOC, T_\text{ref})$,
- $C_1(SOC, T_\text{ref})$,

all defined at a reference temperature $T_\text{ref} = 20^\circ\mathrm{C}$.

The **temperature dependence** of $R_0$, $R_1$, $C_1$ is modeled via an Arrhenius law:

$$
k(T) = k_{\text{ref}} \exp\left(
-\frac{E_a}{R_g}
\left(
\frac{1}{T_K} - \frac{1}{T_{\text{ref},K}}
\right)
\right),
$$

where:

- $k \in \{ R_0, R_1, C_1 \}$,
- $k_{\text{ref}}$ is taken from the $20^\circ\mathrm{C}$ table at the given $SOC$,
- $E_a$ is the activation energy,
- $R_g$ is the gas constant,
- $T_K = T + 273.15$ and $T_{\text{ref},K} = T_\text{ref} + 273.15$.

The ECM’s terminal voltage model is:

$$
V = OCV(SOC) + H_\text{dir}(SOC) + V_{RC} - I R_0(SOC, T)
$$

with $H_\text{dir}(SOC)$ chosen from the charge/discharge hysteresis tables based on current direction.

---

## 2. Process and Measurement Models for the UKF

### 2.1 Process Model

The UKF process model implements:

1. **SOC dynamics**

   Let $Q_\text{cell}$ be the capacity in Ah and $Q_C = 3600 Q_\text{cell}$ be the capacity in Coulombs.

   With $I>0$ = discharge:

   $$ \frac{d(SOC)}{dt} = -\frac{I}{Q_C}.  $$

   Discretized with step $\Delta t$:

   $$ SOC_{k+1} = SOC_k - \frac{\Delta t}{Q_C} I_k.  $$

2. **RC branch dynamics**

   For the diffusion branch:

   $$ \frac{dV_{RC}}{dt} = -\frac{1}{R_1 C_1} V_{RC} + \frac{1}{C_1} I.  $$

   With forward Euler:

   $$ V_{RC,k+1} = V_{RC,k} + \Delta t \left( -\frac{V_{RC,k}}{R_1 C_1} + \frac{I_k}{C_1} \right).  $$

   Here $R_1$ and $C_1$ are looked up from tables and scaled by Arrhenius using the **model ECM** `g_ecm_model`.

3. **Hysteresis dynamics**

   The UKF state $H$ is driven toward a quasi-static hysteresis limit $H_\infty$:

   $$ \frac{dH}{dt} = \frac{H_\infty - H}{\tau_H}, $$

   with $H_\infty$ selected from:

   - $H_\text{dsg}(SOC)$ during discharge ($I>0$),
   - $H_\text{chg}(SOC)$ during charge ($I<0$),
   - or kept near its current value at rest ($I\approx 0$).

   Discrete update:

   $$ H_{k+1} = H_k + \Delta t \, \frac{H_\infty(SOC_k, I_k) - H_k}{\tau_H}.  $$

4. **Thermal dynamics**

   With a lumped thermal model:

   $$ \frac{dT}{dt} = \frac{1}{C_{th}} \left( I^2 R_0(SOC, T) - \frac{T - T_{amb}}{R_{th}} \right), $$

   discretized as:

   $$ T_{k+1} = T_k + \Delta t \frac{ I_k^2 R_0(SOC_k, T_k) - \dfrac{T_k - T_{amb}}{R_{th}} }{C_{th}}.  $$

All required parameters ($R_0$, $R_1$, $C_1$, hysteresis tables, capacity, and thermal parameters) come from `g_ecm_model`.

---

### 2.2 Measurement Model

The measurement model corresponds to the ECM’s terminal voltage:

$$
V_k = OCV(SOC_k)
      + H_k
      + V_{RC,k}
      - I_k R_0(SOC_k, T_k).
$$

The UKF measurement model `ukf_ecm_measurement_model()`:

- Reads $SOC$, $V_{RC}$, $H$, $T$ from the state,
- Uses `ecm_lookup_ocv()` and `ecm_lookup_r0()` from `g_ecm_model`,
- Uses a **global variable** `g_current` for $I_k$ (set prior to each update).

---

## 3. True ECM vs Model ECM

The unit test uses two ECM instances:

1. **`e_true`** — the *plant*:
   - Initialized with `ecm_init_default()`,
   - Tables for $R_0$, $R_1$, $C_1$ are scaled to represent a different, more resistive cell:
     - $R_0$ is multiplied (e.g. by $1.4$),
     - $R_1$ decreased,
     - $C_1$ increased.
   - This ECM is stepped forward with `ecm_step(&e_true, I, T_amb, dt)` and its terminal voltage is computed with `ecm_terminal_voltage(&e_true, I)`.

2. **`g_ecm_model`** — the model parameters used by the UKF:
   - Copy of `ecm_init_default()` (initial guess),
   - Passed to the UKF process and measurement models for parameter lookups,
   - Optionally updated online (if you call `ecm_update_from_measurement()` and `ecm_update_h_from_ukf()`).

In `ukf_ecm_test.c`:

- The **plant** is `e_true`,
- The **UKF dynamics and measurement** use `g_ecm_model` and the UKF state.

---

## 4. Simulation Scenario

The test runs for $N=400$ steps, with time step $\Delta t = 1\ \text{s}$ and ambient temperature $T_{amb}=20^\circ\mathrm{C}$.

The current profile $I_k$ is:

- $k = 0 \ldots 149$: **discharge** at $I = +1.0\ \mathrm{A}$,
- $k = 150 \ldots 199$: **rest** at $I = 0.0\ \mathrm{A}$,
- $k = 200 \ldots 299$: **charge** at $I = -0.5\ \mathrm{A}$,
- $k = 300 \ldots 399$: **rest** at $I = 0.0\ \mathrm{A}$.

At each step:

1. The **true ECM** is stepped:

   ```c
   ecm_step(&e_true, I, T_amb, dt);
   double V_true_clean = ecm_terminal_voltage(&e_true, I);
   ```

2. A noisy measurement is generated:

   $$ V_\text{meas} = V_\text{true} + n_V $$,  $$ n_V \tilde \N(0, \sigma^2_V) $$,

with $\sigma_V \approx 10\ \mathrm{mV}$.

3. The UKF prediction step is executed:

```c
double u[2] = { I, T_amb };
ukf_predict(&ukf, u, dt);
```

4. The UKF update step is executed using $V_\text{meas}$:

```c
g_current = I;   /* for measurement model */
double z[1] = { V_meas };
ukf_update(&ukf, z);
```

5. Optional ECM parameter/hysteresis updates from UKF estimates:

```c
double soc_est = ukf.x[0];
double vrc_est = ukf.x[1];
double H_est   = ukf.x[2];

ecm_update_from_measurement(&g_ecm_model, I, V_meas, vrc_est, dt);

int is_chg = (I < -1e-3) ? 1 : 0;
ecm_update_h_from_ukf(&g_ecm_model, soc_est, H_est, is_chg);

```

6. The UKF’s model-based voltage prediction is computed:

```c
double z_est[1];
ukf_ecm_measurement_model(ukf.x, z_est);
double V_est = z_est[0];
```

7. A CSV line is printed containing:

Time index and current,

* $SOC_\text{true}$ and $SOC_\text{est}$,

* $T_\text{true}$ and $T_\text{est}$,

* $H_\text{est}$,

* $V_\text{true}$, $V_\text{meas}$, $V_\text{est}$,

* True vs model $R_0$, $R_1$, $C_1$ (with Arrhenius scaling).

---

## 5. Output and Interpretation

The CSV header is:
```
k,t,I,
SOC_true,SOC_est,
T_true,T_est,
H_est,
V_true,V_meas,V_est,
R0_true,R0_model,
R1_true,R1_model,
C1_true,C1_model
```

### 5.1 Voltage Tracking

5.1 Voltage Tracking

Plotting $V_\text{true}$, $V_\text{meas}$, and $V_\text{est}$ vs time should show:

* $V_\text{meas}$ oscillating around $V_\text{true}$ with measurement noise,

* $V_\text{est}$ converging close to $V_\text{true}$, indicating:

- The UKF has tuned $SOC$, $V_{RC}$, $H$, and $T$ to be physically consistent,

- The ECM parameter lookups and Arrhenius scaling are being used correctly.


5.2 State Estimation

Plotting $SOC_\text{true}$ vs $SOC_\text{est}$:

* The estimate starts with an initial error,

* During dynamic phases (discharge/charge), the UKF refines $SOC$ based on voltage measurements and the ECM model.

Similarly, $T_\text{true}$ vs $T_\text{est}$ shows how well the lumped thermal dynamics are captured.


5.3 Parameter Behavior

The columns:
* $R0_\text{true}$ vs $R0_\text{model}$,

* $R1_\text{true}$ vs $R1_\text{model}$,

* $C1_\text{true}$ vs $C1_\text{model}$,

indicate whether the model ECM is reasonably consistent with the true cell at the estimated $SOC$ and $T$.

If ecm_update_from_measurement() is enabled, you should see gradual evolution of $R_0$, $R_1$, and $C_1$ in the model towards the true parameters, especially around rest intervals where the identification logic has rich information (voltage steps and exponential decays).

---

## 6. How to Build and Extend

To build:

```bash
gcc -std=c11 -O2 ukf.c ecm.c ukf_ecm_test.c -o ukf_ecm_test -lm
./ukf_ecm_test > ukf_ecm_output.csv
```

ou can then plot with Python, MATLAB, or similar tools.

To extend:

* Replace the simple default tables with real 20-point OCV / hysteresis / impedance data,

* Tune $Q$, $R$, and hysteresis/thermal time constants for your cell,

* Run the same structure with real logged current/voltage/temperature instead of synthetic profiles.

This combined UKF + ECM scaffolding is directly suitable as a BMS fuel-gauge prototype, with clear upgrade paths for more complex models (multiple RC branches, aging parameters, etc.).


