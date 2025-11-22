#!/usr/bin/env python3
"""
Plot the output of ukf_ecm_test.

Usage:

  # 1) Compile and run C test, then plot directly:
  #    gcc -std=c11 -O2 ukf.c ecm.c ukf_ecm_test.c -o ukf_ecm_test -lm
  #    python3 plot_ukf_ecm.py
  #
  #    (This will run ./ukf_ecm_test and plot from its stdout.)

  # 2) Or save output first, then plot from file:
  #    ./ukf_ecm_test > ukf_ecm_output.csv
  #    python3 plot_ukf_ecm.py ukf_ecm_output.csv
"""

import sys
import csv
import subprocess
import matplotlib.pyplot as plt


def get_lines():
    """Return lines either from a file (argv[1]) or from running ./ukf_ecm_test."""
    if len(sys.argv) > 1:
        fname = sys.argv[1]
        with open(fname, "r") as f:
            return f.read().splitlines()
    # No file provided — run the binary
    result = subprocess.run(
        ["./ukf_ecm_test"], check=True, capture_output=True, text=True
    )
    return result.stdout.splitlines()


def parse_csv_lines(lines):
    """Parse CSV lines from ukf_ecm_test into a dict of lists."""
    # Filter out any empty lines
    non_empty = [ln for ln in lines if ln.strip()]
    reader = csv.DictReader(non_empty)

    data = {
        "k": [],
        "t": [],
        "I": [],
        "SOC_true": [],
        "SOC_est": [],
        "T_true": [],
        "T_est": [],
        "H_est": [],
        "V_true": [],
        "V_meas": [],
        "V_est": [],
        "R0_true": [],
        "R0_model": [],
        "R1_true": [],
        "R1_model": [],
        "C1_true": [],
        "C1_model": [],
    }

    for row in reader:
        data["k"].append(int(row["k"]))
        data["t"].append(float(row["t"]))
        data["I"].append(float(row["I"]))

        data["SOC_true"].append(float(row["SOC_true"]))
        data["SOC_est"].append(float(row["SOC_est"]))

        data["T_true"].append(float(row["T_true"]))
        data["T_est"].append(float(row["T_est"]))

        data["H_est"].append(float(row["H_est"]))

        data["V_true"].append(float(row["V_true"]))
        data["V_meas"].append(float(row["V_meas"]))
        data["V_est"].append(float(row["V_est"]))

        data["R0_true"].append(float(row["R0_true"]))
        data["R0_model"].append(float(row["R0_model"]))

        data["R1_true"].append(float(row["R1_true"]))
        data["R1_model"].append(float(row["R1_model"]))

        data["C1_true"].append(float(row["C1_true"]))
        data["C1_model"].append(float(row["C1_model"]))

    return data


def main():
    lines = get_lines()
    data = parse_csv_lines(lines)

    t = data["t"]

    # --- Figure 1: Terminal voltage ---
    fig1, ax1 = plt.subplots(figsize=(9, 4))
    ax1.plot(t, data["V_true"], label="V_true")
    ax1.plot(t, data["V_meas"], label="V_meas", linestyle="none", marker=".", markersize=2)
    ax1.plot(t, data["V_est"], label="V_est (UKF + ECM)", linestyle="--")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Terminal Voltage [V]")
    ax1.set_title("Terminal Voltage: True vs Measured vs Estimated")
    ax1.grid(True)
    ax1.legend()

    # --- Figure 2: SOC and Temperature ---
    fig2, ax2 = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    # SOC
    ax2[0].plot(t, data["SOC_true"], label="SOC_true")
    ax2[0].plot(t, data["SOC_est"], label="SOC_est (UKF)", linestyle="--")
    ax2[0].set_ylabel("SOC [-]")
    ax2[0].set_title("State of Charge and Temperature")
    ax2[0].grid(True)
    ax2[0].legend()

    # Temperature
    ax2[1].plot(t, data["T_true"], label="T_true")
    ax2[1].plot(t, data["T_est"], label="T_est (UKF)", linestyle="--")
    ax2[1].set_xlabel("Time [s]")
    ax2[1].set_ylabel("Temperature [°C]")
    ax2[1].grid(True)
    ax2[1].legend()

    # --- Figure 3: R0, R1, C1 true vs model ---
    fig3, ax3 = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    # R0
    ax3[0].plot(t, data["R0_true"], label="R0_true")
    ax3[0].plot(t, data["R0_model"], label="R0_model", linestyle="--")
    ax3[0].set_ylabel("R0 [Ω]")
    ax3[0].set_title("ECM Parameters: True vs Model")
    ax3[0].grid(True)
    ax3[0].legend()

    # R1
    ax3[1].plot(t, data["R1_true"], label="R1_true")
    ax3[1].plot(t, data["R1_model"], label="R1_model", linestyle="--")
    ax3[1].set_ylabel("R1 [Ω]")
    ax3[1].grid(True)
    ax3[1].legend()

    # C1
    ax3[2].plot(t, data["C1_true"], label="C1_true")
    ax3[2].plot(t, data["C1_model"], label="C1_model", linestyle="--")
    ax3[2].set_xlabel("Time [s]")
    ax3[2].set_ylabel("C1 [F]")
    ax3[2].grid(True)
    ax3[2].legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

