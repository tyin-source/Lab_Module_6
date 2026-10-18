# Phys 39 Module 6, Part 6: integral windup with a low PWM limit (group G)
#
# The one-lump PI run of python/module6_part6_ki_sweep.py (Kp = 8, Ki = 0.36) with the
# PWM limit lowered to u_max = 30, so the initial request Kp * 10.16 = 81 PWM saturates.
# Run twice: with conditional integration (anti-windup on) and with ordinary
# integration (off), as with the "Prevent integral windup" box in the simulation.
#
# Usage (from the repository root):
#   python python/module6_part6_windup.py
# Figure: docs/figures/module_06/part6_G_windup.png

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_part6_ki_sweep import BAND, KP, T_SET, T_AMB, CHI_H, one_lump_pi

REPO_ROOT = Path(__file__).resolve().parents[1]
FIG = REPO_ROOT / "docs" / "figures" / "module_06" / "part6_G_windup.png"
KI = 0.36          # PWM/(C s)
U_MAX = 30.0       # PWM
T_END = 400.0      # s


def main():
    u_ss = (T_SET - T_AMB) / CHI_H
    print(f"Kp {KP:g}, Ki {KI:g}, u_max {U_MAX:g} PWM; u_ss = {u_ss:.2f} PWM, "
          f"q_ss = u_ss/Ki = {u_ss / KI:.0f} C s")
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True, sharey="row",
                             gridspec_kw={"height_ratios": [1.3, 1]})
    for col, aw in enumerate((True, False)):
        t, T, u, u_i = one_lump_pi(KI, u_max=U_MAX, anti_windup=aw, t_end=T_END)
        u_p = KP * (T_SET - T)
        at_limit = abs(u) > U_MAX - 0.01      # conditional integration holds u just at the limit
        t_sat = t[np.flatnonzero(at_limit)[-1]]
        k_pk, k_ui = int(np.argmax(T)), int(np.argmax(u_i))
        outside = np.flatnonzero(abs(T - T_SET) > BAND)
        t_settle = t[outside[-1] + 1]
        k_x = int(np.argmax(T >= T_SET))
        k_q = int(np.argmax(u_i > 0))
        name = "anti-windup on" if aw else "anti-windup off"
        print(f"{name}:")
        print(f"  q starts growing at {t[k_q]:.1f} s (T = {T[k_q]:.2f} C); T first reaches "
              f"T_set at {t[k_x]:.0f} s; u drops below the limit at {t_sat:.0f} s")
        print(f"  u_I peaks at {u_i[k_ui]:.1f} PWM at {t[k_ui]:.0f} s (q = {u_i[k_ui] / KI:.0f} C s)")
        print(f"  T_max {T[k_pk]:.2f} C at {t[k_pk]:.0f} s (overshoot {T[k_pk] - T_SET:+.2f} C); "
              f"settles within +/-{BAND:g} C at {t_settle:.0f} s")

        ax, axu = axes[0, col], axes[1, col]
        ax.plot(t, T, "tab:blue", lw=1.6, label="$T$")
        ax.axhline(T_SET, color="k", ls="--", lw=1, label=f"$T_{{set}}$ = {T_SET:g} °C")
        ax.plot(t[k_pk], T[k_pk], "o", color="tab:blue", ms=6)
        ax.annotate(f"{T[k_pk]:.2f} °C at {t[k_pk]:.0f} s", (t[k_pk], T[k_pk]),
                    xytext=(8, 6) if aw else (10, -12), textcoords="offset points",
                    fontsize=9, color="tab:blue")
        ax.set_title(f"{name}: settles (±{BAND:g} °C) at {t_settle:.0f} s", fontsize=10)
        ax.grid(alpha=0.3)
        axu.plot(t, u_p, "tab:blue", lw=1.2, label="$u_P$")
        axu.plot(t, u_i, "tab:orange", lw=1.4, label="$u_I$")
        axu.plot(t, u, "k", lw=1.6, label="applied $u$")
        axu.axhline(U_MAX, color="tab:red", ls=":", lw=1, label=f"$u_{{\\max}}$ = {U_MAX:g}")
        axu.axhline(u_ss, color="gray", ls="--", lw=0.8, label=f"$u_{{ss}}$ = {u_ss:.2f}")
        axu.axhline(0, color="gray", lw=0.6)
        axu.set_xlabel("Model time (s)")
        axu.grid(alpha=0.3)
    axes[0, 0].set_ylabel("Temperature (°C)")
    axes[1, 0].set_ylabel("Signed PWM")
    axes[0, 0].legend(fontsize=8, loc="lower right")
    axes[1, 1].legend(fontsize=8, loc="upper right")
    axes[1, 0].set_xlim(0, T_END)
    fig.suptitle(f"One-lump PI with a PWM limit: $K_p$ = {KP:g}, $K_i$ = {KI:g}, "
                 f"$u_{{\\max}}$ = {U_MAX:g} PWM", fontsize=11)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG}")


if __name__ == "__main__":
    main()
