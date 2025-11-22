#!/usr/bin/env python3
"""
Plot the output of ecm_test.

Usage:
    # Compile first:
    #   gcc -std=c11 -O2 ecm.c ecm_test.c -o ecm_test -lm
    #
    # Option 1: run and plot directly
    #   python3 plot_ecm.py
    #
    # Option 2: save output then plot
    #   ./ecm_test > ecm_output.csv
    #   python3 plot_ecm.py ecm_output.csv
"""

import sys
import csv
import subprocess
import matplotlib.pyplot as plt


def get_lines():
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r") as f:
            return f.read().splitlines()
    result = subprocess.run(
        ["./ecm_test"], check=True, capture_output=True, text=True
    )
    return result.stdout.splitlines()


def parse(lines):
    reader = csv.DictReader(lines)
    data = {
        "t": [],
        "phase": [],
        "I": [],
        "soc_true": [],
        "soc_model": [],
        "V_true": [],
        "V_model": [],
        "R0_true": [],
        "R0_model": [],
        "R1_true": [],
        "R1_model": [],
        "C1_true": [],
        "C1_model": [],
    }
    for row in reader:
        data["phase"].append(int(row["phase"]))
        data["t"].append(float(row["t"]))
        data["I"].append(float(row["I"]))
        data["soc_true"].append(float(row["soc_true"]))
        data["soc_model"].append(float(row["soc_model"]))
        data["V_true"].append(float(row["V_true"]))
        data["V_model"].append(float(row["V_model"]))
        data["R0_true"].append(float(row["R0_true"]))
        data["R0_model"].append(float(row["R0_model"]))
        data["R1_true"].append(float(row["R1_true"]))
        data["R1_model"].append(float(row["R1_model"]))
        data["C1_true"].append(float(row["C1_true"]))
        data["C1_model"].append(float(row["C1_model"]))
    return data


def main():
    lines = get_lines()
    data = parse(lines)

    t = data["t"]

    fig1, ax1 = plt.subplots(figsize=(9, 4))
    ax1.plot(t, data["V_true"], label="V_true")
    ax1.plot(t, data["V_model"], label="V_model", linestyle="--")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Terminal Voltage [V]")
    ax1.set_title("ECM: True vs Model Voltage")
    ax1.grid(True)
    ax1.legend()

    fig2, ax2 = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    ax2[0].plot(t, data["R0_true"], label="R0_true")
    ax2[0].plot(t, data["R0_model"], label="R0_model", linestyle="--")
    ax2[0].set_ylabel("R0 [ohm]")
    ax2[0].grid(True)
    ax2[0].legend()

    ax2[1].plot(t, data["R1_true"], label="R1_true")
    ax2[1].plot(t, data["R1_model"], label="R1_model", linestyle="--")
    ax2[1].set_ylabel("R1 [ohm]")
    ax2[1].grid(True)
    ax2[1].legend()

    ax2[2].plot(t, data["C1_true"], label="C1_true")
    ax2[2].plot(t, data["C1_model"], label="C1_model", linestyle="--")
    ax2[2].set_xlabel("Time [s]")
    ax2[2].set_ylabel("C1 [F]")
    ax2[2].grid(True)
    ax2[2].legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

