# Phys 39 Module 6, Part 7: one experimental Goldilocks run vs the one-lump model
#
# Reads a run file written by module6_pi_control_gui.py, splits the controlled part at
# the setpoint switches, and reports for each step:
#   T_start   mean of the last 30 s before the step (the baseline mean for the first one)
#   T_ss      mean of the last 60 s of the step, and the steady error e_ss = Tset - T_ss
#   overshoot how far T went past T_ss (above for an upward step, below for a downward one)
#   t_settle  time after the step until T stays within +/-0.1 C of T_ss
#   e(60 s)   error 60 s after the step (shows how much the integral still has to remove)
#   u range   smallest and largest signed PWM sent (cooling negative)
#
# The model is the one-lump block, dT/dt = (T_amb + chi u - T) / tau, with tau = 70 s
# (Part 3), chi_h or chi_c from the Module 4 data, and T_amb = the mean of the baseline
# readings before control started. It is driven by the GUI's own pi_control() with the
# run's gains, setpoint switches, and measured dt; the command is the same integer PWM
# and acts during the next reading, as on the Arduino. It starts from the first
# controlled reading and then runs on its own. Nothing is fitted.
#
# Usage (from the repository root):
#   python python/module6_part7_runs.py data/module_06/<run>.csv --run 1
# Figure: docs/figures/module_06/part7_run<N>.png

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_common import MODULE6_FIG_DIR, TAU_S, module4_susceptibility, read_run
from module6_pi_control_gui import pi_control

BEFORE_S = 30.0     # readings before a step that give T_start
FINAL_S = 60.0      # end of a step that gives T_ss
BAND_C = 0.1        # settling band around T_ss


def one_lump_model(run, t_amb, chi_h, chi_c, substeps=20):
    """Model temperature and commands at the controlled readings of a run."""
    on = run["control_on"] == 1
    tset, dt = run["setpoint_C"][on], run["dt_s"][on]
    kp, ki = run["Kp"][on][0], run["Ki"][on][0]
    u_max, aw = run["u_max"][on][0], bool(run["anti_windup"][on][0])
    n = len(tset)
    T = np.empty(n)
    u = np.empty(n)
    u_p = np.empty(n)
    u_i = np.empty(n)
    T[0], q = run["temperature_C"][on][0], 0.0
    for k in range(n):
        step = pi_control(tset[k], T[k], q, dt[k], kp, ki, u_max, aw)
        q = step.q
        u[k] = step.pwm if step.direction == "HEAT" else -step.pwm
        u_p[k], u_i[k] = step.u_p, step.u_i
        if k + 1 < n:                      # this command acts until the next reading
            chi = chi_h if u[k] >= 0 else chi_c
            T_next = T[k]
            h = dt[k + 1] / substeps
            for _ in range(substeps):
                T_next += h / TAU_S * (t_amb + chi * u[k] - T_next)
            T[k + 1] = T_next
    return T, u, u_p, u_i


def step_metrics(t, T, t_step, t_end, T_start, rising):
    seg = (t >= t_step) & (t <= t_end)
    ts, Ts = t[seg] - t_step, T[seg]
    T_ss = float(Ts[ts >= ts[-1] - FINAL_S].mean())
    over = (Ts.max() - T_ss) if rising else (T_ss - Ts.min())
    outside = np.flatnonzero(abs(Ts - T_ss) > BAND_C)
    if len(outside) == 0:
        t_settle = 0.0
    elif outside[-1] + 1 < len(ts):
        t_settle = float(ts[outside[-1] + 1])
    else:                                   # still outside the band when the step ends
        t_settle = float("nan")
    T_60 = float(np.interp(60.0, ts, Ts))
    return T_ss, max(float(over), 0.0), t_settle, float(ts[-1]), T_60


def main():
    ap = argparse.ArgumentParser(description="Part 7 run vs the one-lump model")
    ap.add_argument("run_file", type=Path)
    ap.add_argument("--run", type=int, required=True, help="Goldilocks run number (1-6)")
    args = ap.parse_args()

    run = read_run(args.run_file)
    chi_h, chi_c = module4_susceptibility()
    on = run["control_on"] == 1
    t0 = run["time_s"][on][0]
    t_all = run["time_s"] - t0
    t, T, tset = t_all[on], run["temperature_C"][on], run["setpoint_C"][on]
    kp, ki = run["Kp"][on][0], run["Ki"][on][0]
    u_max, aw = run["u_max"][on][0], int(run["anti_windup"][on][0])
    sign = np.where(run["dir_cmd"][on] == "HEAT", 1.0, -1.0)
    u_sent = sign * run["pwm_cmd"][on]
    t_amb = float(run["temperature_C"][~on].mean())
    Tm, um, upm, uim = one_lump_model(run, t_amb, chi_h, chi_c)

    print(f"{args.run_file.name}: run {args.run}, Kp {kp:g}, Ki {ki:g}, u_max {u_max:g}, "
          f"anti-windup {'on' if aw else 'off'}; T_amb (baseline) {t_amb:.2f} C, "
          f"{len(t)} controlled readings, dt {run['dt_s'][on][1:].mean():.3f} s")
    switches = np.flatnonzero(np.diff(tset) != 0) + 1
    edges = [0, *switches, len(t)]
    print(" step              Tset  T_start  T_ss meas/model  e_ss meas/model  overshoot meas/model"
          "  settle(+/-0.1) meas/model  e(60 s) meas/model  held   u sent")
    rows = []
    for j in range(len(edges) - 1):
        a, b = edges[j], edges[j + 1] - 1
        t_step, t_end = t[a], t[b]
        if j == 0:
            T_start = t_amb
            name = "from room"
            rising = tset[a] > T_start
        else:
            before = (t >= t_step - BEFORE_S) & (t < t_step)
            T_start = float(T[before].mean())
            rising = tset[a] > tset[a - 1]     # direction of the setpoint change
            name = "up" if rising else "down"
        meas = step_metrics(t, T, t_step, t_end, T_start, rising)
        model = step_metrics(t, Tm, t_step, t_end, T_start, rising)
        seg = slice(a, b + 1)
        rows.append((b, tset[a], meas))
        print(f" {j + 1} {name:9s} {tset[a]:7.2f} {T_start:7.2f}   {meas[0]:6.2f} /{model[0]:6.2f}"
              f"   {tset[a] - meas[0]:+6.2f} /{tset[a] - model[0]:+6.2f}"
              f"     {meas[1]:5.2f} /{model[1]:5.2f}       {meas[2]:5.0f} /{model[2]:5.0f} s"
              f"        {tset[a] - meas[4]:+5.2f} /{tset[a] - model[4]:+5.2f}"
              f"     {meas[3]:4.0f} s  {u_sent[seg].min():+4.0f}..{u_sent[seg].max():+4.0f}")
    L = kp * chi_h
    for level in sorted(set(tset)) if ki == 0 else ():
        print(f"  algebraic one-lump P level for Tset {level:g}: "
              f"T_amb + L/(1+L) (Tset - T_amb) = {t_amb + L / (1 + L) * (level - t_amb):.2f} C "
              f"(L = Kp chi_h = {L:.2f})")

    kind = "P" if ki == 0 else "PI"
    fig, (ax, axu, axe) = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.4, 1]})
    base = ~on & (t_all >= -BEFORE_S)
    ax.plot(t_all[base], run["temperature_C"][base], ".", color="gray", ms=2.5,
            label="baseline (control off)")
    ax.plot(t, T, "k.", ms=2.5, label="measured")
    ax.plot(t, Tm, color="tab:red", lw=1.5,
            label=f"one-lump model (τ = {TAU_S:g} s, χ$_h$ = {chi_h:.3f}, χ$_c$ = {chi_c:.3f}, no fit)")
    ax.step(t, tset, "g--", where="post", lw=1.2, label="$T_{set}$")
    ax.axhline(t_amb, color="gray", ls=":", lw=1, label=f"$T_{{amb}}$ = {t_amb:.2f} °C (baseline)")
    for b, level, meas in rows:           # steady error at the end of each step
        ax.annotate(f"$e_{{ss}}$ = {level - meas[0]:+.2f} °C", (t[b], meas[0]),
                    xytext=(-70, -14), textcoords="offset points", fontsize=8)
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"Run {args.run}: {kind}, $K_p$ = {kp:g} PWM/°C, $K_i$ = {ki:g} PWM/(°C s), "
                 f"$u_{{max}}$ = {u_max:g}, anti-windup {'on' if aw else 'off'}\n{args.run_file.name}",
                 fontsize=10)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.3)
    axu.step(t, u_sent, "k", where="post", lw=1.2, label="PWM sent (signed)")
    axu.plot(t, run["u_P"][on], color="tab:blue", lw=1, label="$u_P$")
    if ki > 0:
        axu.plot(t, run["u_I"][on], color="tab:orange", lw=1.4, label="$u_I$")
    axu.step(t, um, color="tab:red", where="post", lw=1, alpha=0.8, label="model PWM")
    axu.axhline(0, color="gray", lw=0.6)
    axu.set_ylabel("Signed PWM")
    axu.legend(fontsize=7.5, loc="upper right", ncol=2)
    axu.grid(alpha=0.3)
    axe.plot(t, run["error_C"][on], color="m", lw=1.2, label="measured")
    axe.plot(t, tset - Tm, color="tab:red", lw=1, label="model")
    axe.axhline(0, color="gray", lw=0.6)
    axe.set_ylabel("$e = T_{set} - T$ (°C)")
    axe.set_xlabel("Time since control started (s)")
    axe.legend(fontsize=7.5, loc="upper right")
    axe.grid(alpha=0.3)
    fig.tight_layout()
    MODULE6_FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = MODULE6_FIG_DIR / f"part7_run{args.run}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
