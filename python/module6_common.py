# Phys 39 Module 6 - shared helpers
# (copied from module5_common.py; Module 6 adds the time constant tau and the
# one-lump PI damping ratio, and reads the Part 7 run files)
#
# Used by the Part 7 GUI and analysis scripts:
#   - Module 4 open-loop susceptibilities (from data/module_04/part4_steady_state.csv)
#   - the Part 3 time constant tau and the one-lump PI damping ratio
#   - reading one run file written by module6_pi_control_gui.py
#
# Susceptibility convention (Module 6 handout): u is the signed PWM command,
# u > 0 heats and u < 0 cools, and T_ss = T_amb + chi u in both directions with
# the susceptibility of the active direction, so chi_h and chi_c are both
# positive. (Module 5 wrote the cooling value as |chi_c|; the number is the same.)
# The Module 4 fit used signed PWM (cooling PWM negative), so its cooling
# slope m_c is already the positive chi_c.

import csv
import math
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE4_DATA = REPO_ROOT / "data" / "module_04" / "part4_steady_state.csv"
MODULE6_DATA_DIR = REPO_ROOT / "data" / "module_06"
MODULE6_FIG_DIR = REPO_ROOT / "docs" / "figures" / "module_06"

TAU_S = 70.0      # one-lump time constant (s), Module 6 Part 3 open-loop steps (70 +/- 3 s)


def module4_susceptibility(path: Path = MODULE4_DATA):
    """Return (chi_h, chi_c) in C per PWM count from the Module 4 Part 4 data.

    Same straight-line fits as Module 4's part4_temp_vs_pwm.py: heating uses the
    PWM 0 point and the HEAT points; cooling uses the COOL points and PWM 0.
    """
    signed_pwm = []
    temp = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            signed_pwm.append(float(row["signed_pwm"]))
            temp.append(float(row["steady_temperature_C"]))
    signed_pwm = np.array(signed_pwm)
    temp = np.array(temp)

    heat = signed_pwm >= 0
    cool = signed_pwm <= 0
    chi_h, _ = np.polyfit(signed_pwm[heat], temp[heat], 1)
    chi_c, _ = np.polyfit(signed_pwm[cool], temp[cool], 1)
    return float(chi_h), float(chi_c)


def directional_susceptibility(t_set: float, t_amb: float):
    """Pick the Module 4 susceptibility for the direction the controller must
    drive to reach t_set from t_amb.

    Returns (direction, chi) with direction "HEAT" or "COOL".
    """
    chi_h, chi_c = module4_susceptibility()
    if t_set >= t_amb:
        return "HEAT", chi_h
    return "COOL", chi_c


def one_lump_pi_zeta(chi: float, kp: float, ki: float, tau: float = TAU_S):
    """One-lump PI damping ratio zeta = (1 + chi Kp) / (2 sqrt(tau chi Ki)).

    Returns inf for Ki = 0 (P only: a single real pole, no oscillation).
    The step response does not overshoot when Ki <= Kp / tau (Module 6 Part 6).
    """
    if ki <= 0.0:
        return math.inf
    return (1.0 + chi * kp) / (2.0 * math.sqrt(tau * chi * ki))


def read_run(path: Path):
    """Read one run CSV from module6_pi_control_gui.py into numpy arrays."""
    cols = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for name in reader.fieldnames:
            cols[name] = []
        for row in reader:
            for name in reader.fieldnames:
                cols[name].append(row[name])

    text_cols = {"dir_cmd"}
    out = {}
    for name, values in cols.items():
        if name in text_cols:
            out[name] = np.array(values)
        else:
            out[name] = np.array([float(v) for v in values])
    return out
