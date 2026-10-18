# Phys 39 Module 6, Part 7: the six experimental Goldilocks runs side by side
#
# Overlays the matched 30 -> 35 C (upward) and 35 -> 30 C (downward) steps of the six
# retained runs, each plotted against the time since its own setpoint switch, and
# prints the step metrics of module6_part7_runs.py for all of them.
#
# Usage (from the repository root):
#   python python/module6_part7_summary.py
# Figure: docs/figures/module_06/part7_goldilocks_summary.png

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_common import MODULE6_DATA_DIR, MODULE6_FIG_DIR, read_run
from module6_part7_runs import BEFORE_S, step_metrics

RUNS = [  # (run, file, label, colour)
    (1, "p_20261018-045404_Kp2_Ki0_T30-35.csv", "P, $K_p$ = 2", "#9ecae1"),
    (2, "p_20261018-054301_Kp8_Ki0_T30-35.csv", "P, $K_p$ = 8", "#3182bd"),
    (3, "p_20261018-062146_Kp32_Ki0_T30-35.csv", "P, $K_p$ = 32", "#08306b"),
    (4, "pi_20261018-063831_Kp6_Ki0.04_T30-35.csv", "PI, $K_p$ = 6, $K_i$ = 0.04", "#fdae6b"),
    (5, "pi_20261018-072600_Kp6_Ki0.08_T30-35.csv", "PI, $K_p$ = 6, $K_i$ = 0.08", "#e6550d"),
    (6, "pi_20261018-083656_Kp6_Ki0.5_T30-35.csv", "PI, $K_p$ = 6, $K_i$ = 0.5", "#7f2704"),
]
SHOW_S = 400.0      # seconds after each step shown


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    print(" run  step   Tset  e_ss (C)  overshoot (C)  settle +/-0.1 C (s)  e(60 s) (C)")
    for run, name, label, colour in RUNS:
        r = read_run(MODULE6_DATA_DIR / name)
        on = r["control_on"] == 1
        t = r["time_s"][on] - r["time_s"][on][0]
        T, tset = r["temperature_C"][on], r["setpoint_C"][on]
        switches = np.flatnonzero(np.diff(tset) != 0) + 1          # up, then down
        edges = [*switches, len(t)]
        for ax, j, step in ((axes[0], 0, "up"), (axes[1], 1, "down")):
            a, b = edges[j], edges[j + 1] - 1
            T_start = float(T[(t >= t[a] - BEFORE_S) & (t < t[a])].mean())
            T_ss, over, t_settle, _, T_60 = step_metrics(t, T, t[a], t[b], T_start,
                                                       rising=(step == "up"))
            print(f"  {run}   {step:4s}  {tset[a]:5.1f}   {tset[a] - T_ss:+6.2f}       "
                  f"{over:5.2f}          {t_settle:5.0f}            {tset[a] - T_60:+5.2f}")
            sel = (t >= t[a]) & (t <= min(t[b], t[a] + SHOW_S))
            ax.plot(t[sel] - t[a], T[sel], color=colour, lw=1.6, label=f"Run {run}: {label}")
    for ax, title, level in ((axes[0], "Upward step, 30 → 35 °C", 35.0),
                             (axes[1], "Downward step, 35 → 30 °C", 30.0)):
        ax.axhline(level, color="g", ls="--", lw=1, label="$T_{set}$")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Time since the setpoint switch (s)")
        ax.set_xlim(0, SHOW_S)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Temperature (°C)")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Part 7: six experimental Goldilocks runs (one apparatus, $u_{max}$ = 255, "
                 "anti-windup on)", fontsize=11)
    fig.tight_layout()
    out = MODULE6_FIG_DIR / "part7_goldilocks_summary.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
