# Phys 39 Module 6, Part 5: the Module 5 high-gain P trace vs the one-lump P model
#
# Marks the late-time temperature T_ss and the first overshoot of the measured
# Kp = 32 run and overlays the one-lump prediction with no fitted parameter:
#   T_{n+1} = T_n + (dt / tau) * (T_amb + chi * u_n - T_n),   u_n = clamp(Kp (T_set - T_n))
# with chi_h = 0.485 C/PWM (chi_c = 0.172 if u < 0), tau = 70 s (Part 3), T_amb = the
# baseline mean before P control started, and T(0) = the first reading under control.
#
# Usage (from the repository root):
#   python archive/python/module6_part5_highgain.py
# Figure: docs/figures/module_06/part5_highgain_Kp32_vs_one_lump.png

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
RUN = REPO_ROOT / "data" / "module_05" / "p_20261014-120316_Tset30.0_Kp32.csv"
FIG = REPO_ROOT / "docs" / "figures" / "module_06" / "part5_highgain_Kp32_vs_one_lump.png"
CHI_H, CHI_C = 0.485, 0.172    # C/PWM, Module 4
TAU = 70.0                     # s, Part 3
PWM_MAX = 255.0
SS_WINDOW_S = 120.0            # late-time window for T_ss (as in Module 5)


def read_run(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    col = lambda k: np.array([float(r[k]) for r in rows])
    heat = col("heat_cool") == 1
    return (col("time_s"), col("temperature_C"), col("setpoint_C"), col("Kp"),
            np.where(heat, col("pwm_arduino"), -col("pwm_arduino")),
            col("control_on").astype(bool))


def one_lump_p(t_end, T0, t_amb, t_set, kp, dt=0.01):
    n = int(round(t_end / dt))
    T = np.empty(n + 1)
    u = np.empty(n + 1)
    T[0] = T0
    for i in range(n + 1):
        u[i] = min(PWM_MAX, max(-PWM_MAX, kp * (t_set - T[i])))
        if i < n:
            chi = CHI_H if u[i] >= 0 else CHI_C
            T[i + 1] = T[i] + dt / TAU * (t_amb + chi * u[i] - T[i])
    return np.arange(n + 1) * dt, T, u


def main():
    t, T, t_set_col, kp_col, u, on = read_run(RUN)
    t0 = t[on][0]
    rel = t - t0
    t_set, kp = float(t_set_col[on][0]), float(kp_col[on][0])
    t_amb = float(T[~on].mean())
    tc, Tc = rel[on], T[on]
    T_ss = float(Tc[tc >= tc[-1] - SS_WINDOW_S].mean())
    L = CHI_H * kp
    T_ss_model = t_amb + L * (t_set - t_amb) / (1 + L)

    k_cross = int(np.argmax(Tc >= T_ss))
    k_peak = int(np.argmax(Tc))
    after = (tc > tc[k_peak]) & (tc < tc[k_peak] + 25)
    k_dip = int(np.flatnonzero(after)[np.argmin(Tc[after])])
    tm, Tm, um = one_lump_p(tc[-1], float(Tc[0]), t_amb, t_set, kp)

    print(f"{RUN.name}: Kp {kp:g}, T_set {t_set:g} C, baseline T_amb {t_amb:.2f} C, "
          f"T(0) {Tc[0]:.2f} C")
    print(f"  late-time T_ss {T_ss:.2f} C (last {SS_WINDOW_S:g} s); one-lump T_ss "
          f"{T_ss_model:.2f} C, tau_cl {TAU / (1 + L):.2f} s")
    print(f"  first crossing of T_ss at {tc[k_cross]:.1f} s; first overshoot {Tc[k_peak]:.2f} C "
          f"at {tc[k_peak]:.1f} s (+{Tc[k_peak] - T_ss:.2f} C above T_ss, "
          f"{100 * (Tc[k_peak] - T_ss) / (T_ss - Tc[0]):.0f} % of the step; "
          f"{Tc[k_peak] - t_set:+.2f} C past T_set)")
    print(f"  following dip {Tc[k_dip]:.2f} C at {tc[k_dip]:.1f} s ({Tc[k_dip] - T_ss:+.2f} C)")
    print(f"  one-lump model: max {Tm.max():.3f} C (never above its T_ss {T_ss_model:.3f} C)")

    win = (rel >= -5) & (rel <= 80)
    fig, (ax, axp) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    ax.plot(rel[win], T[win], "k.-", lw=1, ms=3.5, label="measured, Module 5 (14 Oct 2026)")
    sel = tm <= 80
    ax.plot(tm[sel], Tm[sel], color="tab:red", lw=1.8,
            label=f"one-lump P model: τ = {TAU:g} s, χ$_h$ = {CHI_H}, no fit")
    ax.axhline(t_set, color="tab:green", ls="--", lw=1, label=f"$T_{{set}}$ = {t_set:g} °C")
    ax.axhline(T_ss, color="k", ls=":", lw=1.2, label=f"late-time $T_{{ss}}$ = {T_ss:.2f} °C")
    ax.plot(tc[k_cross], Tc[k_cross], "o", mfc="none", mec="tab:blue", ms=9)
    ax.annotate(f"crosses $T_{{ss}}$\nat {tc[k_cross]:.0f} s", (tc[k_cross], Tc[k_cross]),
                xytext=(-75, 25), textcoords="offset points", color="tab:blue", fontsize=9,
                arrowprops=dict(arrowstyle="->", color="tab:blue"))
    ax.plot(tc[k_peak], Tc[k_peak], "o", color="tab:blue", ms=7)
    ax.annotate(f"first overshoot {Tc[k_peak]:.2f} °C at {tc[k_peak]:.0f} s\n"
                f"+{Tc[k_peak] - T_ss:.2f} °C above $T_{{ss}}$",
                (tc[k_peak], Tc[k_peak]), xytext=(18, -22), textcoords="offset points",
                color="tab:blue", fontsize=9)
    ax.annotate(f"{Tc[k_dip] - T_ss:+.2f} °C", (tc[k_dip], Tc[k_dip]), xytext=(5, -15),
                textcoords="offset points", color="tab:blue", fontsize=9)
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"$K_p$ = {kp:g} PWM/°C: measured transient vs one-lump P model\n{RUN.name}",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    axp.step(rel[win], u[win], "k", where="pre", lw=1.2,      # each reading's PWM held during its average
             label="measured (applied by the Arduino)")
    axp.plot(tm[sel], um[sel], color="tab:red", lw=1.4, label="one-lump P model")
    axp.axhline(0, color="gray", lw=0.6)
    axp.set_ylabel("Signed PWM")
    axp.set_xlabel("Time since P control started (s)")
    axp.legend(fontsize=8)
    axp.grid(alpha=0.3)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG}")


if __name__ == "__main__":
    main()
