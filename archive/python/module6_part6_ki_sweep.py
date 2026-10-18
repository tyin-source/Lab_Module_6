# Phys 39 Module 6, Part 6: one-lump PI step response for a sweep of Ki (group F)
#
# Reproduces the course simulation (one_lump, measured, pi, anti-windup on) with the
# measured-parameter form of the model, so C, H, and P_u never appear:
#   q_{n+1} = q_n + e_n dt          (skipped while it would drive u further into the limit)
#   u_n     = clamp(Kp e_n + Ki q_{n+1}),   e_n = T_set - T_n
#   T_{n+1} = T_n + (dt / tau) * (T_amb + chi * u_n - T_n)
# with chi = chi_h for u >= 0 and chi_c for u < 0. The closed loop is
#   dT/dt = (T_amb + chi (Kp e + Ki q) - T) / tau,   dq/dt = e,
# a second-order system with damping ratio zeta = (1 + chi Kp) / (2 sqrt(tau chi Ki))
# and a zero at s = -Ki/Kp from the PI controller.
#
# Usage (from the repository root):
#   python archive/python/module6_part6_ki_sweep.py
# Figures: docs/figures/module_06/part6_F_ki_sweep.png and part6_F_pole_zero.png (Ki = Kp/tau,
# where the PI zero cancels the slow pole and the step response is first order)

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG = REPO_ROOT / "docs" / "figures" / "module_06" / "part6_F_ki_sweep.png"
FIG_PZ = REPO_ROOT / "docs" / "figures" / "module_06" / "part6_F_pole_zero.png"
CHI_H, CHI_C = 0.485, 0.172    # C/PWM, Module 4
TAU = 70.0                     # s, Part 3
T_AMB = T0 = 19.84             # C, as entered in the simulation
T_SET = 30.0                   # C
KP = 8.0                       # PWM/C
KIS = (0.08, 0.18, 0.28, 0.36, 0.70, 2.0)   # PWM/(C s)
PWM_MAX = 255.0
DT = 0.05                      # s, as in the simulation
T_END = 600.0                  # s
BAND = 0.2                     # C, settling band (2 % of the 10.16 C step)


def one_lump_pi(ki, kp=KP, u_max=PWM_MAX, anti_windup=True, dt=DT, t_end=T_END, T_start=T0):
    """Return time, temperature, applied command, and the I term u_I = Ki q (q starts at 0)."""
    n = int(round(t_end / dt))
    T = np.empty(n + 1)
    u = np.empty(n + 1)
    u_i = np.empty(n + 1)
    T[0], q = T_start, 0.0
    for i in range(n + 1):
        e = T_SET - T[i]
        proposed = kp * e + ki * (q + e * dt)
        if not (anti_windup and ((proposed > u_max and e > 0) or (proposed < -u_max and e < 0))):
            q += e * dt
        u_i[i] = ki * q
        u[i] = min(u_max, max(-u_max, kp * e + u_i[i]))
        if i < n:
            chi = CHI_H if u[i] >= 0 else CHI_C
            T[i + 1] = T[i] + dt / TAU * (T_AMB + chi * u[i] - T[i])
    return np.arange(n + 1) * dt, T, u, u_i


def settling_time(t, T):
    outside = np.flatnonzero(abs(T - T_SET) > BAND)
    return t[outside[-1] + 1] if len(outside) else 0.0


def pole_zero():
    """Ki = Kp/tau puts the PI zero -Ki/Kp on the open-loop pole -1/tau."""
    ki_pz = KP / TAU
    ki_crit = (1 + CHI_H * KP) ** 2 / (4 * TAU * CHI_H)          # zeta = 1
    tau_cl = TAU / (CHI_H * KP)
    print(f"Pole-zero cancellation: Ki = Kp/tau = {ki_pz:.4f}; zeta = 1 at Ki = {ki_crit:.4f}; "
          f"first-order tau_cl = tau/(chi Kp) = {tau_cl:.2f} s, settles in "
          f"{tau_cl * np.log((T_SET - T_AMB) / BAND):.1f} s")
    cases = [(f"$K_i$ = $K_p/τ$ = {ki_pz:.3f}", ki_pz, "tab:green"),
             ("$K_i$ = 0.12 (simulation field)", 0.12, "tab:blue"),
             (f"$K_i$ = {ki_crit:.3f} (ζ = 1)", ki_crit, "tab:red")]
    fig, (ax, axz) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [1.4, 1]})
    for label, ki, color in cases:
        t, T, _, _ = one_lump_pi(ki)
        zeta = (1 + CHI_H * KP) / (2 * np.sqrt(TAU * CHI_H * ki))
        print(f"  Ki {ki:.4f}: zeta {zeta:.3f}, T_max {T.max():.3f} C, "
              f"settles within +/-{BAND:g} C at {settling_time(t, T):.1f} s")
        ax.plot(t, T, color=color, lw=1.6, label=f"{label}, ζ = {zeta:.2f}")
        axz.plot(t, T - T_SET, color=color, lw=1.6)
        if ki == ki_pz:
            T_fo = T_AMB + (T_SET - T_AMB) * (1 - np.exp(-t / tau_cl))
            print(f"  Ki = Kp/tau vs first order: largest difference {abs(T - T_fo).max():.4f} C")
    # Exact cancellation the simulation's fields accept (Kp in steps of 1, Ki in steps of 0.02)
    t7, T7, _, _ = one_lump_pi(0.10, kp=7.0)
    tau_cl7 = TAU / (CHI_H * 7.0)
    T7_fo = T_AMB + (T_SET - T_AMB) * (1 - np.exp(-t7 / tau_cl7))
    print(f"  Kp 7, Ki 0.10 (= Kp/tau): zeta {(1 + CHI_H * 7) / (2 * np.sqrt(TAU * CHI_H * 0.10)):.3f}, "
          f"T_max {T7.max():.3f} C, settles at {settling_time(t7, T7):.1f} s; tau_cl {tau_cl7:.2f} s, "
          f"largest difference from first order {abs(T7 - T7_fo).max():.4f} C")
    # The cancelled pole is still in the loop: Zero integral at the PI steady state
    # (T = T_set, q = 0) recovers on the open-loop time scale.
    for ki in (ki_pz, 0.36, 0.08):
        tz, Tz, _, _ = one_lump_pi(ki, T_start=T_SET)
        print(f"  Zero integral at steady state, Ki {ki:.4f}: dip to {Tz.min():.2f} C at "
              f"{tz[np.argmin(Tz)]:.0f} s, back within +/-{BAND:g} C at {settling_time(tz, Tz):.0f} s")
    ax.plot(t, T_fo, "k:", lw=1.6, label=f"first order, $τ_{{cl}}$ = τ/(χ$_h$$K_p$) = {tau_cl:.1f} s")
    ax.axhline(T_SET, color="k", ls="--", lw=0.8)
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"One-lump PI, $K_p$ = {KP:g} PWM/°C: where the overshoot disappears", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    axz.axhspan(-BAND, BAND, color="gray", alpha=0.2, label=f"±{BAND:g} °C")
    axz.axhline(0, color="k", lw=0.8)
    axz.set_ylim(-0.6, 0.6)
    axz.set_xlim(0, 250)
    axz.set_ylabel("$T - T_{set}$ (°C), zoomed")
    axz.set_xlabel("Model time (s)")
    axz.legend(fontsize=8, loc="lower right")
    axz.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_PZ, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG_PZ}")


def main():
    u_ss = (T_SET - T_AMB) / CHI_H
    print(f"Kp {KP:g} PWM/C, T_amb = T0 = {T_AMB} C, T_set {T_SET:g} C, tau {TAU:g} s, "
          f"chi_h {CHI_H}; u_ss = {u_ss:.2f} PWM; settling band +/-{BAND:g} C")
    print(" Ki     zeta   T_max   overshoot  t_peak  t_settle  min u   area to first T_set   u_ss/Ki")
    fig, (ax, axu) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1.2]})
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(KIS)))
    for ki, color in zip(KIS, colors):
        t, T, u, _ = one_lump_pi(ki)
        zeta = (1 + CHI_H * KP) / (2 * np.sqrt(TAU * CHI_H * ki))
        k_pk = int(np.argmax(T))
        t_settle = settling_time(t, T)
        k_x = int(np.argmax(T >= T_SET))         # first reach of T_set (0 if never)
        area = np.sum(T_SET - T[:k_x]) * DT if k_x else np.sum(T_SET - T) * DT
        print(f" {ki:<5g}  {zeta:.2f}  {T[k_pk]:6.2f}  {T[k_pk] - T_SET:+8.2f}  "
              f"{t[k_pk]:6.0f}  {t_settle:7.0f}  {u.min():6.1f}  "
              f"{area:12.0f} C s{'' if k_x else ' (never)':8s}  {u_ss / ki:7.0f} C s")
        ax.plot(t, T, color=color, lw=1.6, label=f"$K_i$ = {ki:g}, ζ = {zeta:.2f}")
        axu.plot(t, u, color=color, lw=1.3)
    ax.axhline(T_SET, color="k", ls="--", lw=1, label=f"$T_{{set}}$ = {T_SET:g} °C")
    ax.axhspan(T_SET - BAND, T_SET + BAND, color="gray", alpha=0.2,
               label=f"±{BAND:g} °C settling band")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"One-lump PI, $K_p$ = {KP:g} PWM/°C, $K_i$ sweep "
                 f"(τ = {TAU:g} s, χ$_h$ = {CHI_H}, χ$_c$ = {CHI_C})", fontsize=10)
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    ax.grid(alpha=0.3)
    axu.axhline(u_ss, color="k", ls=":", lw=1, label=f"$u_{{ss}}$ = {u_ss:.2f} PWM")
    axu.axhline(0, color="gray", lw=0.6)
    axu.set_ylabel("Signed PWM $u$")
    axu.set_xlabel("Model time (s)")
    axu.set_xlim(0, 250)
    axu.legend(fontsize=8)
    axu.grid(alpha=0.3)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG}")
    pole_zero()


if __name__ == "__main__":
    main()
