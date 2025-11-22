# Generic Unscented Kalman Filter (UKF) in C

This repository provides a **fully configurable**, **embedded-friendly**, **Unscented Kalman Filter (UKF)** implementation in C, along with a sample unit test and visualization script.  
It supports nonlinear systems, bounded RAM usage, and runs without external libraries.


## 📌 Overview

The Unscented Kalman Filter (UKF) is an estimation algorithm useful for tracking hidden states in nonlinear systems. Examples include:

- Battery state estimation
- Robotics and SLAM
- Navigation, IMU sensor fusion
- Automotive and industrial control

UKF improves over the Extended Kalman Filter (EKF) by **avoiding Jacobians** and instead using **deterministic sampling (sigma points)** to capture how uncertainty flows through the nonlinear model.


## 📖 Mathematical Model

The UKF assumes a **nonlinear state-space system**:

### State (process) model

\[
\mathbf{x}_{k+1} = f(\mathbf{x}_k, \mathbf{u}_k, w_k)
\]

### Measurement model

\[
\mathbf{z}_k = h(\mathbf{x}_k, v_k)
\]

Where:

| Symbol | Meaning |
|--------|---------|
| \( \mathbf{x}_k \) | System state (dimension \(n_x\)) |
| \( \mathbf{z}_k \) | Measurement vector (dimension \(n_z\)) |
| \( \mathbf{u}_k \) | Control input (optional) |
| \( w_k, v_k \) | Zero-mean process and measurement noise |
| \( Q, R \) | Noise covariances |


## ✨ Sigma Points & Scaling

Rather than linearizing, the UKF approximates mean and covariance by evaluating sample points around the current estimate. These are called **sigma points**.

### Sigma Point Generation

For a state mean \( \mathbf{x} \) and covariance \(P\):

\[
\lambda = \alpha^2 (n + \kappa) - n, \qquad \gamma = \sqrt{n + \lambda}
\]

\[
\chi_0 = \mathbf{x}
\]

\[
\chi_i = \mathbf{x} + \gamma \mathbf{l}_i \quad (i = 1..n)
\]

\[
\chi_{i+n} = \mathbf{x} - \gamma \mathbf{l}_i \quad (i = 1..n)
\]

where \( \mathbf{l}_i \) is the *i-th column* of the Cholesky factor \(L\):

\[
P = L L^T
\]


### Weighting

Mean weights:

\[
W_0^{(m)} = \frac{\lambda}{n + \lambda}, \qquad 
W_i^{(m)} = \frac{1}{2(n + \lambda)}
\]

Covariance weights:

\[
W_0^{(c)} = \frac{\lambda}{n + \lambda} + (1 - \alpha^2 + \beta)
\]
\[
W_i^{(c)} = \frac{1}{2(n + \lambda)}
\]

### Meaning of Parameters

| Parameter | Role | Recommended |
|-----------|------|-------------|
| \( \alpha \) | Spread of sigma points | \(10^{-3}\) to \(10^{-1}\) |
| \( \beta \) | Encodes prior knowledge (Gaussian → 2) | 2 |
| \( \kappa \) | Secondary scaling | 0 or 3-n |


## 🔮 Prediction Step

Given state \( \mathbf{x}_k \), covariance \(P_k\), and input \( \mathbf{u}_k \):

1. Generate sigma points from \( (\mathbf{x}_k, P_k) \)
2. Propagate through process model:

\[
\chi_{k+1|k}^i = f(\chi_k^i, \mathbf{u}_k, \Delta t)
\]

3. Predicted mean:

\[
\hat{\mathbf{x}}_{k+1|k} = \sum W_i^{(m)} \chi_{k+1|k}^i
\]

4. Predicted covariance:

\[
P_{k+1|k} = Q + \sum W_i^{(c)} (\chi_{k+1|k}^i - \hat{\mathbf{x}}_{k+1|k})(\chi_{k+1|k}^i - \hat{\mathbf{x}}_{k+1|k})^T
\]


## 🛰️ Update Step

Given a measurement \( \mathbf{z}_{k+1} \):

1. Transform sigma points into measurement space:

\[
\mathbf{z}_{k+1}^i = h(\chi_{k+1|k}^i)
\]

2. Predicted measurement:

\[
\hat{\mathbf{z}}_{k+1} = \sum W_i^{(m)} \mathbf{z}_{k+1}^i
\]

3. Measurement covariance:

\[
S = R + \sum W_i^{(c)} (\mathbf{z}_{k+1}^i - \hat{\mathbf{z}}_{k+1})(\mathbf{z}_{k+1}^i - \hat{\mathbf{z}}_{k+1})^T
\]

4. Cross-covariance:

\[
P_{xz} = \sum W_i^{(c)} (\chi_{k+1|k}^i - \hat{\mathbf{x}}_{k+1|k})(\mathbf{z}_{k+1}^i - \hat{\mathbf{z}}_{k+1})^T
\]

5. Kalman gain:

\[
K = P_{xz} S^{-1}
\]

6. State update:

\[
\mathbf{x}_{k+1} = \hat{\mathbf{x}}_{k+1|k} + K(\mathbf{z}_{k+1} - \hat{\mathbf{z}}_{k+1})
\]

7. Covariance update:

\[
P_{k+1} = P_{k+1|k} - K S K^T
\]


## 🛠️ Usage (C API Summary)

```c
ukf_init(&ukf, n_x, n_z, alpha, beta, kappa);
ukf_set_models(&ukf, my_fx, my_hx);
ukf_set_state(&ukf, x0, P0);
ukf_set_noise(&ukf, Q, R);

ukf_predict(&ukf, u_k, dt);

if (have_measurement)
    ukf_update(&ukf, z_k);
````

## 📂 Contents

| File               | Description                                  |
| ------------------ | -------------------------------------------- |
| `ukf.h`            | Public API                                   |
| `ukf.c`            | Implementation (no malloc, no external deps) |
| `ukf_test.c`       | Example: constant-velocity tracking          |
| `plot_ukf_test.py` | Python plotting of test output               |


## 🧪 Example Result

After running `ukf_test`, you’ll see:

* Normal convergence behavior
* Handling of missing measurements
* Behavior when measurement noise is (almost) zero

The corresponding plots are generated using:

```bash
python3 plot_ukf_test.py
```


## 🏷️ License

MIT License — free to use in research, commercial products, and embedded firmware.

```

