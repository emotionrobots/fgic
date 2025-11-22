# Generic Unscented Kalman Filter (UKF) in C

This repository provides a **fully configurable**, **embedded-friendly**, **Unscented Kalman Filter (UKF)** implementation in C, along with:

- A detailed theory of operation,
- A reference unit test,
- A Python plotting utility.

The implementation is designed for **microcontrollers and real-time systems**:  
🟢 **No dynamic allocation**  
🟢 **No external dependencies**  
🟢 **Small RAM footprint**

---

## 📌 Overview

The Unscented Kalman Filter (UKF) is used to estimate the internal states of **nonlinear systems**. Common applications include:

- Battery management (SOC, SOH, hysteresis, temperature)
- Robotics and navigation (SLAM, IMU fusion)
- Signal processing
- Automotive and industrial control

Unlike the Extended Kalman Filter (EKF), the UKF does **not compute Jacobians**.  
Instead, it uses **deterministic sigma points** to propagate uncertainty.

---

## 📖 Mathematical Model

The UKF assumes a nonlinear state-space system:

### 🔸 State (process) model

$$
\mathbf{x}_{k+1} = f(\mathbf{x}_k, \mathbf{u}_k, w_k)
$$

### 🔹 Measurement model

$$
\mathbf{z}_k = h(\mathbf{x}_k, v_k)
$$

Where:

| Symbol | Meaning |
|--------|---------|
| $\mathbf{x}_k$ | System state (size $n_x$) |
| $\mathbf{z}_k$ | Measurement (size $n_z$) |
| $\mathbf{u}_k$ | Known input |
| $w_k, v_k$ | Process / measurement noise |
| $Q, R$ | Noise covariances |

---

## ✨ Sigma Points & Scaling

To approximate nonlinear uncertainty, the filter generates **$2n + 1$ deterministic points**:

$$
\lambda = \alpha^2 (n + \kappa) - n
$$

$$
\gamma = \sqrt{n + \lambda}
$$

$$
\chi_0 = \mathbf{x}
$$

$$
\chi_i = \mathbf{x} + \gamma \mathbf{l}_i \quad (i = 1..n)
$$

$$
\chi_{i+n} = \mathbf{x} - \gamma \mathbf{l}_i \quad (i = 1..n)
$$

Where $\mathbf{l}_i$ are columns of the Cholesky factor $L$:

$$
P = L L^T
$$

---

### Weighting

**Mean weights**

$$
W_0^{(m)} = \frac{\lambda}{n + \lambda}, \qquad 
W_i^{(m)} = \frac{1}{2(n + \lambda)}
$$

**Covariance weights**

$$
W_0^{(c)} = \frac{\lambda}{n + \lambda} + (1 - \alpha^2 + \beta)
$$

$$
W_i^{(c)} = \frac{1}{2(n + \lambda)}
$$

---

### Parameter Guidelines

| Parameter | Purpose | Suggested |
|-----------|----------|----------|
| $\alpha$ | Sigma spread | $10^{-3}$ to $10^{-1}$ |
| $\beta$ | Distribution knowledge | $2$ for Gaussian |
| $\kappa$ | Secondary tuning | $0$ or $3-n$ |

---

## 🔮 Prediction Step

1. Generate sigma points  
2. Propagate through nonlinear model $f$  
3. Compute predicted mean and covariance

$$
\hat{\mathbf{x}}_{k+1|k} = \sum W_i^{(m)} f(\chi_k^i)
$$

$$
P_{k+1|k} = Q + \sum W_i^{(c)} (\chi_{k+1|k}^i - \hat{\mathbf{x}}_{k+1|k})(\chi_{k+1|k}^i - \hat{\mathbf{x}}_{k+1|k})^T
$$

---

## 🛰️ Update Step

Given a measurement $\mathbf{z}_{k+1}$:

1. Transform sigma points through measurement model $h$
2. Compute predicted measurement and covariance
3. Compute cross-covariance
4. Apply Kalman gain

$$
K = P_{xz} S^{-1}
$$

5. Update state and covariance

$$
\mathbf{x}_{k+1} = \hat{\mathbf{x}}_{k+1|k} + K(\mathbf{z}_{k+1} - \hat{\mathbf{z}}_{k+1})
$$

$$
P_{k+1} = P_{k+1|k} - K S K^T
$$

---

## 🧪 The UKF Unit Test Explained

`ukf_test.c` demonstrates filter behavior using a **1D constant-velocity model**:

### 🚗 System Model

- State:

$$
\mathbf{x} = 
\begin{bmatrix}
\text{position} \\
\text{velocity}
\end{bmatrix}
$$

- Process model:

$$
x_{k+1} = x_k + v_k \Delta t, \qquad
v_{k+1} = v_k
$$

- Measurement model:

$$
z_k = \text{position}
$$

---

### 📌 **Case 1 — Normal Operation**

> Measurement is available every update.

✔ Estimate converges from a wrong initial guess  
✔ Velocity estimate improves **despite never being measured**  
✔ Covariance decreases steadily

---

### 🚫 **Case 2 — Missing Measurements**

> From step 20 to 29, no measurement is provided.

Behavior:

- Covariance grows
- Position and velocity drift slowly
- When measurements resume, the filter **corrects instantly**

---

### 🎯 **Case 3 — Very Small Measurement Noise**

> Measurement noise covariance $R \approx 0$

Behavior:

- State estimate hugs the measurement
- Dynamics contribute very little
- Shows how $R$ controls **measurement dominance**

---

## 📊 Plotting the Output

Compile and visualize:

```bash
gcc -std=c11 -O2 ukf.c ukf_test.c -o ukf_test -lm
python3 plot_ukf_test.py
```

# 🏷️  License

MIT License — free to use in research, commercial products, and embedded firmware.

