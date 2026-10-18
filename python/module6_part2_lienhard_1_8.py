# Phys 39 Module 6, Part II pre-class: Lienhard Problem 1.8 (lumped cooling of a copper sphere)
#
# A copper sphere, D = 2.5 cm, starts at 40 C in a 0 C air stream with h = 15 W/(m^2 K).
# Lumped energy balance:  rho c V dT/dt = -h A (T - T_inf)
#   =>  T(t) = T_inf + (T_0 - T_inf) exp(-t / tau),   tau = rho c V / (h A) = rho c D / (6 h)
# Biot number Bi = h (V/A) / k = h D / (6 k).
# Copper properties near room temperature (Lienhard, Appendix A): rho = 8954 kg/m^3,
# c = 384 J/(kg K), k = 398 W/(m K).
#
# Usage (from the repository root):
#   python python/module6_part2_lienhard_1_8.py
# Figure: docs/figures/module_06/part2_lienhard_1_8.png

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
FIG = REPO_ROOT / "docs" / "figures" / "module_06" / "part2_lienhard_1_8.png"
RHO, C_P, K = 8954.0, 384.0, 398.0      # copper: kg/m^3, J/(kg K), W/(m K)
D, H_CONV = 0.025, 15.0                 # m, W/(m^2 K)
T0, T_INF = 40.0, 0.0                   # C


def main():
    V_over_A = D / 6
    tau = RHO * C_P * V_over_A / H_CONV
    bi = H_CONV * V_over_A / K
    C = RHO * C_P * np.pi * D**3 / 6      # heat capacity, J/K
    hA = H_CONV * np.pi * D**2            # conductance to the air, W/K
    print(f"V/A = {V_over_A * 1e3:.3f} mm, Bi = {bi:.2e}, C = {C:.1f} J/K, hA = {hA:.4f} W/K, "
          f"tau = C/(hA) = {tau:.0f} s = {tau / 60:.1f} min")
    for frac in (0.5, 0.1, 0.01):
        print(f"  T - T_inf falls to {frac:g} of its start after {tau * np.log(1 / frac):.0f} s")
    print(f"  heat released 40 -> 0 C: {C * (T0 - T_INF):.0f} J")

    t = np.linspace(0, 5 * tau, 500)
    T = T_INF + (T0 - T_INF) * np.exp(-t / tau)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t / 60, T, "tab:red", lw=2, label=f"$T(t) = 40\\,e^{{-t/\\tau}}$ °C, τ = {tau:.0f} s")
    ax.axhline(T0 * np.exp(-1), color="gray", ls=":", lw=1)
    ax.axvline(tau / 60, color="gray", ls=":", lw=1)
    ax.annotate(f"t = τ: {T0 * np.exp(-1):.1f} °C", (tau / 60, T0 * np.exp(-1)), xytext=(8, 8),
                textcoords="offset points", fontsize=9)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Sphere temperature (°C)")
    ax.set_title(f"Lienhard Problem 1.8: lumped cooling of a copper sphere (Bi = {bi:.1e})",
                 fontsize=10)
    ax.set_xlim(0, 5 * tau / 60)
    ax.set_ylim(0, 41)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=120)
    plt.close(fig)
    print(f"Saved {FIG}")


if __name__ == "__main__":
    main()
