#!/usr/bin/env python3
"""
Plot the output of ukf_test.c

Usage:
  1) Compile and run directly:
       gcc -std=c11 -O2 ukf.c ukf_test.c -o ukf_test -lm
       python3 plot_ukf_test.py          # runs ./ukf_test and plots
  2) Or save the output to a file:
       ./ukf_test > ukf_output.txt
       python3 plot_ukf_test.py ukf_output.txt
"""

import sys
import re
import subprocess
import matplotlib.pyplot as plt


def get_ukf_output():
    """Get ukf_test output either from a file (if provided) or by running ./ukf_test."""
    if len(sys.argv) > 1:
        fname = sys.argv[1]
        with open(fname, "r") as f:
            return f.read().splitlines()
    # Run ./ukf_test and capture stdout
    result = subprocess.run(
        ["./ukf_test"], check=True, capture_output=True, text=True
    )
    return result.stdout.splitlines()


def parse_output(lines):
    """
    Parse ukf_test.c output into three cases.

    Case 1 lines:
      k, t, x_true, v_true, z_meas, x_est, v_est

    Case 2 lines:
      k=.. t=.. (with meas|NO meas) x_true=.. x_est=..

    Case 3 lines:
      k=.. t=.. x_true=.. z_meas=.. x_est=..
    """

    # Storage
    case1 = {"t": [], "x_true": [], "v_true": [], "z_meas": [], "x_est": [], "v_est": []}
    case2 = {"t": [], "x_true": [], "x_est": [], "has_meas": []}
    case3 = {"t": [], "x_true": [], "z_meas": [], "x_est": []}

    current_case = None

    # Regexes for Cases 2 and 3
    re_case2 = re.compile(
        r"k=\s*(\d+)\s+t=([0-9.eE+-]+)\s+\((with meas|NO meas)\)\s+"
        r"x_true=([0-9.eE+-]+)\s+x_est=([0-9.eE+-]+)"
    )
    re_case3 = re.compile(
        r"k=\s*(\d+)\s+t=([0-9.eE+-]+)\s+"
        r"x_true=([0-9.eE+-]+)\s+z_meas=([0-9.eE+-]+)\s+x_est=([0-9.eE+-]+)"
    )

    for line in lines:
        line_stripped = line.strip()

        if not line_stripped:
            continue

        if "Case 1:" in line_stripped:
            current_case = 1
            continue
        if "Case 2:" in line_stripped:
            current_case = 2
            continue
        if "Case 3:" in line_stripped:
            current_case = 3
            continue

        # Skip header line in Case 1
        if current_case == 1 and line_stripped.startswith("k,"):
            continue

        if current_case == 1:
            # Example: " 0, 0.10, 0.1000, 1.0000, 0.1243, 0.0200, 0.5000"
            if "," in line_stripped:
                parts = [p.strip() for p in line_stripped.split(",")]
                if len(parts) == 7:
                    # k = int(parts[0])  # not used
                    t = float(parts[1])
                    x_true = float(parts[2])
                    v_true = float(parts[3])
                    z_meas = float(parts[4])
                    x_est = float(parts[5])
                    v_est = float(parts[6])

                    case1["t"].append(t)
                    case1["x_true"].append(x_true)
                    case1["v_true"].append(v_true)
                    case1["z_meas"].append(z_meas)
                    case1["x_est"].append(x_est)
                    case1["v_est"].append(v_est)

        elif current_case == 2:
            m = re_case2.search(line_stripped)
            if m:
                # k = int(m.group(1))  # not used
                t = float(m.group(2))
                meas_flag = m.group(3)
                x_true = float(m.group(4))
                x_est = float(m.group(5))
                has_meas = "with meas" in meas_flag

                case2["t"].append(t)
                case2["x_true"].append(x_true)
                case2["x_est"].append(x_est)
                case2["has_meas"].append(has_meas)

        elif current_case == 3:
            m = re_case3.search(line_stripped)
            if m:
                # k = int(m.group(1))  # not used
                t = float(m.group(2))
                x_true = float(m.group(3))
                z_meas = float(m.group(4))
                x_est = float(m.group(5))

                case3["t"].append(t)
                case3["x_true"].append(x_true)
                case3["z_meas"].append(z_meas)
                case3["x_est"].append(x_est)

    return case1, case2, case3


def plot_case1(case1):
    if not case1["t"]:
        print("Case 1: no data parsed, skipping plot.")
        return

    t = case1["t"]
    x_true = case1["x_true"]
    x_est = case1["x_est"]
    z_meas = case1["z_meas"]
    v_true = case1["v_true"]
    v_est = case1["v_est"]

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    fig.suptitle("Case 1: Normal operation")

    # Position
    ax[0].plot(t, x_true, label="x_true")
    ax[0].plot(t, x_est, label="x_est (UKF)")
    ax[0].scatter(t, z_meas, s=10, marker="x", label="z_meas")
    ax[0].set_ylabel("Position")
    ax[0].grid(True)
    ax[0].legend()

    # Velocity
    ax[1].plot(t, v_true, label="v_true")
    ax[1].plot(t, v_est, label="v_est (UKF)")
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Velocity")
    ax[1].grid(True)
    ax[1].legend()


def plot_case2(case2):
    if not case2["t"]:
        print("Case 2: no data parsed, skipping plot.")
        return

    t = case2["t"]
    x_true = case2["x_true"]
    x_est = case2["x_est"]
    has_meas = case2["has_meas"]

    fig, ax = plt.subplots(figsize=(8, 4))
    fig.suptitle("Case 2: Missing measurements (k=20..29)")

    ax.plot(t, x_true, label="x_true")
    ax.plot(t, x_est, label="x_est (UKF)")

    # Mark points with and without measurement
    t_meas = [ti for ti, hm in zip(t, has_meas) if hm]
    x_meas = [xi for xi, hm in zip(x_est, has_meas) if hm]
    t_nomeas = [ti for ti, hm in zip(t, has_meas) if not hm]
    x_nomeas = [xi for xi, hm in zip(x_est, has_meas) if not hm]

    ax.scatter(t_meas, x_meas, s=25, marker="o", label="update w/ meas")
    ax.scatter(t_nomeas, x_nomeas, s=25, marker="x", label="predict only")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position")
    ax.grid(True)
    ax.legend()


def plot_case3(case3):
    if not case3["t"]:
        print("Case 3: no data parsed, skipping plot.")
        return

    t = case3["t"]
    x_true = case3["x_true"]
    x_est = case3["x_est"]
    z_meas = case3["z_meas"]

    fig, ax = plt.subplots(figsize=(8, 4))
    fig.suptitle("Case 3: Very small measurement noise")

    ax.plot(t, x_true, label="x_true")
    ax.plot(t, x_est, label="x_est (UKF)")
    ax.scatter(t, z_meas, s=25, marker="x", label="z_meas")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position")
    ax.grid(True)
    ax.legend()


def main():
    lines = get_ukf_output()
    case1, case2, case3 = parse_output(lines)

    plot_case1(case1)
    plot_case2(case2)
    plot_case3(case3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

