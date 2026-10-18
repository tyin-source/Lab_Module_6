# Phys 39 Module 6, Part 4: one-lump open-loop model vs the measured open-loop steps
#
# The model is the measured-parameter Euler update of the handout (the same equation the
# course simulation integrates in one_lump / measured / open_loop mode):
#   T_{n+1} = T_n + (dt / tau) * (T_amb + chi * u - T_n)
# with chi = chi_h for u > 0 and chi_c for u < 0. Nothing is fitted to the trace: chi
# comes from Module 4, tau from Part 3, and T_amb is the measured level before the step.
#
# Usage (from the repository root):
#   python python/module6_part4_open_loop.py
# Figure: docs/figures/module_06/part4_open_loop_model_vs_measured.png

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_tau63 import read_trace

REPO_ROOT = Path(__file__).resolve().parents[1]
TRACE = REPO_ROOT / "data" / "module_06" / "module4_run_20261017_164351.csv"
FIG = REPO_ROOT / "docs" / "figures" / "module_06" / "part4_open_loop_model_vs_measured.png"
CHI_H = 0.485          # C/PWM, Module 4
CHI_C = 0.172          # C/PWM, Module 4
TAU = 70.0             # s, Part 3
BASELINE_S = 30.0      # readings before the step that give T_amb


def one_lump(t_end, t_amb, u, dt):
    """Euler integration from T = T_amb at t = 0 with a constant command u."""
    chi = CHI_H if u > 0 else CHI_C
    n = int(round(t_end / dt))
    T = np.empty(n + 1)
    T[0] = t_amb
    for i in range(n):
        T[i + 1] = T[i] + dt / TAU * (t_amb + chi * u - T[i])
    return np.arange(n + 1) * dt, T


def main():
    t, T, u = read_trace(TRACE)
    changes = np.where(np.diff(u) != 0)[0] + 1
    steps = [("Heating step, u = +25 PWM", changes[0], changes[1]),
             ("Cooling step, u = −25 PWM", changes[2], changes[3])]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1]})
    print(f"chi_h {CHI_H}, chi_c {CHI_C} C/PWM, tau {TAU:g} s; trace {TRACE.name}")
    for col, (title, k0, k1) in enumerate(steps):
        cmd = float(u[k0])
        before = (t < t[k0]) & (t >= t[k0] - BASELINE_S)
        t_amb = float(T[before].mean())
        ts = t[k0:k1] - t[k0]
        Ts = T[k0:k1]
        T_ss_model = t_amb + (CHI_H if cmd > 0 else CHI_C) * cmd
        T_ss_meas = float(Ts[ts >= ts[-1] - 60].mean())

        # Delta t check: the model at the reading times for three step sizes
        curves = {}
        for dt in (1.0, 0.1, 0.01):
            tm, Tm = one_lump(ts[-1] + 1, t_amb, cmd, dt)
            curves[dt] = np.interp(ts, tm, Tm)
        model = curves[0.01]
        res = Ts - model
        print(f"{title}: T_amb {t_amb:.2f} C, held {ts[-1]:.0f} s")
        print(f"  T_ss: model {T_ss_model:.2f} C, measured {T_ss_meas:.2f} C "
              f"(difference {T_ss_meas - T_ss_model:+.2f} C; measured step "
              f"{100 * ((T_ss_meas - t_amb) / (T_ss_model - t_amb) - 1):+.1f} % vs model)")
        print(f"  measured - model: rms {res.std():.2f} C, largest {res[np.argmax(abs(res))]:+.2f} C "
              f"at {ts[np.argmax(abs(res))]:.0f} s")
        print(f"  dt check: max |T(dt=1 s) - T(dt=0.01 s)| = {abs(curves[1.0] - model).max():.3f} C, "
              f"max |T(dt=0.1 s) - T(dt=0.01 s)| = {abs(curves[0.1] - model).max():.4f} C")

        ax, axr = axes[0, col], axes[1, col]
        ax.plot(ts, Ts, "k.", ms=2.5, label="measured (17 Oct 2026)")
        ax.plot(ts, model, "tab:red", lw=1.6,
                label=f"one-lump model, χ = {CHI_H if cmd > 0 else CHI_C} °C/PWM, τ = {TAU:g} s")
        ax.axhline(T_ss_model, color="tab:red", ls=":", lw=1,
                   label=f"model $T_{{ss}}$ = {T_ss_model:.2f} °C")
        ax.axhline(T_ss_meas, color="k", ls=":", lw=1,
                   label=f"measured $T_{{ss}}$ = {T_ss_meas:.2f} °C")
        ax.set_title(f"{title} from $T_{{amb}}$ = {t_amb:.2f} °C", fontsize=10)
        ax.set_ylabel("Temperature (°C)")
        ax.legend(fontsize=8, loc="center right")
        ax.grid(alpha=0.3)
        axr.plot(ts, res, "k", lw=1)
        axr.axhline(0, color="tab:red", lw=1)
        axr.set_ylabel("measured − model (°C)")
        axr.set_xlabel("Time since the step (s)")
        axr.grid(alpha=0.3)
    fig.suptitle("Part 4: open-loop one-lump model (no fitted parameters) vs measurement",
                 fontsize=11)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG}")


if __name__ == "__main__":
    main()
