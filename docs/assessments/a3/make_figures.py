# Phys 39 Module 6 -- A3 figures (vector PDF, sized for the REVTeX reprint layout)
#
#   figures/exp_grid.pdf       six experimental Goldilocks runs (Question 2)
#   figures/sim_grid.pdf       six two-lump simulations, same axes (Question 3)
#   figures/model_compare.pdf  Run 3 and Run 6 upward steps: measured, one-lump, two-lump
#                              (Questions 3 and 5)
#
# The experiment comes from data/module_06/ (the six runs listed in
# python/module6_part7_summary.py). The simulations are the ones of
# python/module6_part2_goldilocks.py (course model, G = 2 W/K, C_T/C = 0.25); the
# one-lump model is the one of python/module6_part7_runs.py.
#
# Usage (from the repository root):
#   python docs/assessments/a3/make_figures.py

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_common import MODULE6_DATA_DIR, module4_susceptibility, read_run
from module6_part7_runs import one_lump_model
from module6_part7_summary import RUNS
import module6_part2_goldilocks as P

FIG_DIR = HERE / "figures"
TEXT_W, COL_W = 7.0, 3.4          # inches
BEFORE_S = 30.0
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
                     "lines.linewidth": 1.0, "axes.linewidth": 0.6, "pdf.fonttype": 42})
EXP, SIM, ONE, SET = "k", "tab:blue", "tab:red", "tab:green"


def load(n, name):
    """Measured run and its two-lump simulation, times from the start of control."""
    r = read_run(MODULE6_DATA_DIR / name)
    on = r["control_on"] == 1
    t_all = r["time_s"] - r["time_s"][on][0]
    t, T, tset = t_all[on], r["temperature_C"][on], r["setpoint_C"][on]
    kp, ki = float(r["Kp"][on][0]), float(r["Ki"][on][0])
    t_amb = float(r["temperature_C"][~on].mean())
    starts = [0, *(np.flatnonzero(np.diff(tset) != 0) + 1)]
    c = P.shared_config(t_amb, float(T[0]), kp, ki, float(r["u_max"][on][0]),
                        bool(r["anti_windup"][on][0]))
    s = P.simulate_run(c, "p" if ki == 0 else "pi", t[starts], tset[starts],
                       t[-1] + r["dt_s"][on][-1])
    base = ~on & (t_all >= -BEFORE_S)
    sign = np.where(r["dir_cmd"][on] == "HEAT", 1.0, -1.0)
    return dict(n=n, run=r, c=c, kp=kp, ki=ki, t_amb=t_amb, t=t, T=T, tset=tset,
                t_base=t_all[base], T_base=r["temperature_C"][base],
                u=sign * r["pwm_cmd"][on], u_I=r["u_I"][on], sim=s)


def label(d, sim):
    kind = f"P, $K_p$ = {d['kp']:g}" if d["ki"] == 0 else \
        f"PI, $K_p$ = {d['kp']:g}, $K_i$ = {d['ki']:g}"
    if not sim:
        return f"Run {d['n']}: {kind}"
    if d["ki"] == 0:
        return f"Run {d['n']} sim.: {kind}; $\\zeta_{{2L,P}}$ = {P.M.two_lump_p_damping_ratio(d['c']):.2f}"
    kc = P.M.two_lump_pi_critical_integral_gain(d["c"])
    return (f"Run {d['n']} sim.: {kind}; $K_{{i,crit}}$ = {kc:.2f}, "
            f"$\\rho_I$ = {d['ki'] / kc:.3f}")


def limits(d):
    s = d["sim"]
    temps = np.concatenate([d["T"], d["T_base"], s["Tm"], [d["t_amb"]], d["tset"]])
    pwm = np.concatenate([d["u"], s["u"], [0.0]])
    pad_T, pad_u = 0.05 * np.ptp(temps), 0.08 * np.ptp(pwm)
    return ((d["t_base"][0], d["t"][-1]), (temps.min() - pad_T, temps.max() + pad_T),
            (pwm.min() - pad_u, pwm.max() + pad_u))


def panel(fig, cell, d, sim, letter):
    sub = cell.subgridspec(2, 1, height_ratios=[2, 1], hspace=0.08)
    ax, axu = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    xl, yl, ul = limits(d)
    s = d["sim"]
    if sim:
        ax.plot(d["t"], d["T"], ".", color="0.75", ms=1.2, label="experiment")
        ax.plot(s["t"], s["Tm"], color=SIM, lw=1.1, label="simulation $T_m$")
        axu.step(d["t"], d["u"], color="0.75", where="post", lw=0.8)
        axu.plot(s["t"], s["u"], color=SIM, lw=0.9, label="applied $u$")
        if d["ki"] > 0:
            axu.plot(s["t"], s["u_I"], color="tab:orange", lw=0.8, label="$u_I$")
    else:
        ax.plot(d["t_base"], d["T_base"], ".", color="0.6", ms=1.2)
        ax.plot(d["t"], d["T"], ".", color=EXP, ms=1.2, label="measured $T$")
        axu.step(d["t"], d["u"], color=EXP, where="post", lw=0.8, label="signed PWM sent")
        if d["ki"] > 0:
            axu.plot(d["t"], d["u_I"], color="tab:orange", lw=0.8, label="$u_I$")
    ax.step(d["t"], d["tset"], color=SET, ls="--", where="post", lw=0.9, label="$T_{set}$")
    ax.axhline(d["t_amb"], color="0.5", ls=":", lw=0.7, label="$T_{amb}$")
    axu.axhline(0, color="0.5", lw=0.5)
    for a in (ax, axu):
        a.set_xlim(*xl)
        a.grid(alpha=0.25, lw=0.4)
    ax.set_ylim(*yl)
    axu.set_ylim(*ul)
    ax.tick_params(labelbottom=False)
    ax.set_title(f"({letter}) {label(d, sim)}", loc="left", pad=2)
    ax.set_ylabel("$T$ (°C)")
    axu.set_ylabel("PWM")
    return ax, axu


def grid(runs, sim, out):
    fig = plt.figure(figsize=(TEXT_W, 7.9))
    gs = fig.add_gridspec(3, 2, hspace=0.32, wspace=0.2, left=0.07, right=0.99,
                          top=0.955, bottom=0.05)
    letters = "abcdef"
    handles = {}
    for i, d in enumerate(runs):              # left column P (runs 1-3), right PI (4-6)
        ax, axu = panel(fig, gs[i % 3, i // 3], d, sim, letters[i])
        for a in (ax, axu):
            for h, l in zip(*a.get_legend_handles_labels()):
                handles.setdefault(l, h)
        if i % 3 == 2:
            axu.set_xlabel("Time since control started (s)")
    fig.legend(handles.values(), handles.keys(), loc="upper center", ncol=len(handles),
               frameon=False, bbox_to_anchor=(0.5, 1.0), markerscale=4)
    fig.savefig(FIG_DIR / out)
    plt.close(fig)
    print(f"Saved {FIG_DIR / out}")


def model_compare(runs):
    chi_h, chi_c = module4_susceptibility()
    fig, axes = plt.subplots(2, 1, figsize=(COL_W, 3.5))
    for ax, d, show in ((axes[0], runs[2], 60.0), (axes[1], runs[5], 220.0)):
        T1, _, _, _ = one_lump_model(d["run"], d["t_amb"], chi_h, chi_c)
        t_up = d["t"][np.flatnonzero(np.diff(d["tset"]) > 0)[0] + 1]
        s = d["sim"]
        m = (d["t"] >= t_up - 5) & (d["t"] <= t_up + show)
        ms = (s["t"] >= t_up - 5) & (s["t"] <= t_up + show)
        ax.plot(d["t"][m] - t_up, d["T"][m], ".", color=EXP, ms=2, label="measured")
        ax.plot(d["t"][m] - t_up, T1[m], color=ONE, lw=1, label="one-lump model")
        ax.plot(s["t"][ms] - t_up, s["Tm"][ms], color=SIM, lw=1, label="two-lump $T_m$")
        ax.axhline(35, color=SET, ls="--", lw=0.8)
        ax.set_xlim(-5, show)
        ax.grid(alpha=0.25, lw=0.4)
        ax.set_ylabel("$T$ (°C)")
        ax.set_title(f"({'ab'[d['n'] == 6]}) {label(d, False)}, upward step 30 → 35 °C",
                     loc="left", pad=2)
    axes[0].legend(loc="lower right", markerscale=2)
    axes[1].set_xlabel("Time since the setpoint switch (s)")
    fig.tight_layout(pad=0.3, h_pad=0.6)
    fig.savefig(FIG_DIR / "model_compare.pdf")
    plt.close(fig)
    print(f"Saved {FIG_DIR / 'model_compare.pdf'}")


def main():
    FIG_DIR.mkdir(exist_ok=True)
    runs = [load(n, name) for n, name, _, _ in RUNS]
    grid(runs, False, "exp_grid.pdf")
    grid(runs, True, "sim_grid.pdf")
    model_compare(runs)


if __name__ == "__main__":
    main()
