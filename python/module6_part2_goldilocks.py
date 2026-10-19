# Phys 39 Module 6, Part II: the six two-lump simulation Goldilocks plots
#
# Runs the course two-lump model (Lab_6_7_modeling_tec.py, its own advance_model) for
# the six experimental Goldilocks runs of Part 7 with ONE shared physical parameter set:
#   measured mode: chi_c = 0.172, chi_h = 0.485 C/PWM, tau = 70 s, V_TEC = 10.95 V,
#                  R_M = 1.50 ohm  ->  C, H, P_u,c, r from constants_from_measurements()
#   locked in Part II Step 2 from Run 3 (K_p = 32): G = 2 W/K, C_T/C = 0.25
# Only the controller gains change. Each simulation takes from its run file:
#   T_amb   mean of the baseline readings before control started
#   T_0     the first controlled reading (both lumps start there, as after Reset)
#   the setpoint schedule (the setpoint is switched at the measured times, without
#   resetting the integral), K_p, K_i, u_max, and anti-windup
# The model controller acts on T_m every Euler step (dt = 0.05 s); it has no sampling,
# averaging, or one-line command delay. Nothing is fitted here.
#
# For each step it prints the same metrics for the measurement and the simulation:
#   e_ss      Tset - (mean of the last 60 s of the step)
#   rise      10 % -> 90 % of the change from T_start to T_ss
#   overshoot past T_ss, and the time of that extreme after the step
#   settle    time until T stays within +/-0.1 C of T_ss ("never" if it is still
#             outside the band when the step ends)
# T_start is T_amb for the step from room temperature and the mean of the last 30 s
# before the switch for the others (module6_part7_runs.py).
#
# Usage (from the repository root):
#   python python/module6_part2_goldilocks.py            # all six runs
#   python python/module6_part2_goldilocks.py --run 3    # one run
# Figures: docs/figures/module_06/part2_sim_run<N>.png

import argparse
import sys
import types
from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from module6_common import MODULE6_DATA_DIR, MODULE6_FIG_DIR, read_run
from module6_part7_runs import BEFORE_S, step_metrics
from module6_part7_summary import RUNS

CHI_C, CHI_H, TAU, V_TEC, R_M = 0.172, 0.485, 70.0, 10.95, 1.50   # as typed in the GUI
G_W_PER_C = 2.0          # locked in Part II Step 2
CT_FRACTION = 0.25       # locked in Part II Step 2
DT_S = 0.05


def course_model():
    """Import the course simulation for its model functions (its window is not used)."""
    try:
        import tkinter  # noqa: F401
    except ModuleNotFoundError:          # headless Python: give the GUI imports empty modules
        for name in ("tkinter", "tkinter.ttk", "tkinter.messagebox"):
            sys.modules[name] = types.ModuleType(name)
        sys.modules["tkinter"].Tk = sys.modules["tkinter"].PhotoImage = object
    import Lab_6_7_modeling_tec
    return Lab_6_7_modeling_tec


M = course_model()


def shared_config(t_amb, t0, kp, ki, u_max, anti_windup):
    d = M.constants_from_measurements(CHI_C, CHI_H, TAU, V_TEC, R_M)
    c = replace(M.ModelConfig(), ambient_c=t_amb, initial_c=t0,
                measured_cooling_chi_c_per_pwm=CHI_C, measured_heating_chi_c_per_pwm=CHI_H,
                measured_tau_s=TAU, tec_voltage_v=V_TEC, tec_resistance_ohm=R_M,
                total_capacitance_j_per_c=d.total_capacitance_j_per_c,
                heat_loss_w_per_c=d.heat_loss_w_per_c,
                tec_power_w_per_pwm=d.cooling_power_w_per_pwm,
                heating_to_cooling_ratio=d.heating_to_cooling_ratio,
                coupling_w_per_c=G_W_PER_C, tec_capacitance_fraction=CT_FRACTION,
                kp_pwm_per_c=kp, ki_pwm_per_c_s=ki, pwm_limit=u_max, dt_s=DT_S,
                anti_windup=anti_windup)
    M.validate_config(c)
    return c


def kp_crit(c):
    """Lowest K_p at which the two-lump P loop is underdamped (heating)."""
    c_t, c_m = c.tec_capacitance_j_per_c, c.measured_capacitance_j_per_c
    g, h = c.coupling_w_per_c, c.heat_loss_w_per_c
    return ((c_t * (g + h) + c_m * g) ** 2 / (4 * c_t * c_m * g) - h) / M.heating_tec_coefficient(c)


def simulate_run(c, mode, t_switch, setpoints, t_end):
    """Two-lump model with the setpoint switched at the measured times."""
    state = M.ModelState(0.0, c.initial_c, c.initial_c)
    edges = [*t_switch[1:], t_end]
    out = []
    for sp, stop in zip(setpoints, edges):
        cs = replace(c, setpoint_c=sp)
        while state.time_s < stop - 1e-9:
            r = M.advance_model(state, cs, "two_lump", mode)
            out.append((r.time_s, r.measured_temperature_c, r.tec_temperature_c, sp,
                        r.applied_pwm, r.p_pwm, r.i_pwm))
    return {k: np.array(v) for k, v in zip(("t", "Tm", "T", "sp", "u", "u_P", "u_I"), zip(*out))}


def rise_time(ts, Ts, T_start, T_ss):
    """10 % -> 90 % time of the change from T_start to T_ss (nan if not reached)."""
    frac = (Ts - T_start) / (T_ss - T_start)
    i10, i90 = np.flatnonzero(frac >= 0.1), np.flatnonzero(frac >= 0.9)
    if len(i10) == 0 or len(i90) == 0:
        return float("nan")
    return float(ts[i90[0]] - ts[i10[0]])


def metrics(t, T, t_step, t_end, T_start, rising):
    T_ss, over, t_settle, _, _ = step_metrics(t, T, t_step, t_end, T_start, rising)
    seg = (t >= t_step) & (t <= t_end)
    ts, Ts = t[seg] - t_step, T[seg]
    t_peak = float(ts[np.argmax(Ts) if rising else np.argmin(Ts)])
    return T_ss, rise_time(ts, Ts, T_start, T_ss), over, t_peak, t_settle


def caption(c, kp, ki):
    if ki == 0:
        z = M.two_lump_p_damping_ratio(c)
        side = "above" if kp > kp_crit(c) else "below"
        return (f"two-lump P, $K_p$ = {kp:g} PWM/°C: $\\zeta_{{2L,P}}$ = {z:.3f} "
                f"({M.damping_description(z)}), {side} $K_{{p,crit}}$ = {kp_crit(c):.2f} PWM/°C")
    kic = M.two_lump_pi_critical_integral_gain(c)
    return (f"two-lump PI, $K_p$ = {kp:g} PWM/°C, $K_i$ = {ki:g} PWM/(°C s): "
            f"$K_{{i,crit}}$ = {kic:.2f} PWM/(°C s), $\\rho_I$ = {ki / kic:.3f}")


def do_run(n, name):
    run = read_run(MODULE6_DATA_DIR / name)
    on = run["control_on"] == 1
    t0 = run["time_s"][on][0]
    t_all = run["time_s"] - t0
    t, T, tset = t_all[on], run["temperature_C"][on], run["setpoint_C"][on]
    kp, ki = float(run["Kp"][on][0]), float(run["Ki"][on][0])
    u_max, aw = float(run["u_max"][on][0]), bool(run["anti_windup"][on][0])
    sign = np.where(run["dir_cmd"][on] == "HEAT", 1.0, -1.0)
    u_sent = sign * run["pwm_cmd"][on]
    t_amb = float(run["temperature_C"][~on].mean())
    switches = np.flatnonzero(np.diff(tset) != 0) + 1
    starts = [0, *switches]
    mode = "p" if ki == 0 else "pi"
    c = shared_config(t_amb, float(T[0]), kp, ki, u_max, aw)
    s = simulate_run(c, mode, t[starts], tset[starts], t[-1] + run["dt_s"][on][-1])

    print(f"\nRun {n}: {name}; T_amb {t_amb:.2f} C, T_0 {T[0]:.2f} C; "
          + caption(c, kp, ki).replace("$", "").replace("\\", "").replace("{", "").replace("}", ""))
    print("  step        Tset   e_ss meas/sim   rise meas/sim (s)   overshoot meas/sim (C)"
          "   at meas/sim (s)   settle meas/sim (s)")
    edges = [*starts, len(t)]
    for j in range(len(edges) - 1):
        a, b = edges[j], edges[j + 1] - 1
        if j == 0:
            name_j, rising = "from room", tset[a] > t_amb
            T_start_m = T_start_s = t_amb
        else:
            rising = tset[a] > tset[a - 1]
            name_j = "up" if rising else "down"
            T_start_m = float(T[(t >= t[a] - BEFORE_S) & (t < t[a])].mean())
            T_start_s = float(s["Tm"][(s["t"] >= t[a] - BEFORE_S) & (s["t"] < t[a])].mean())
        m = metrics(t, T, t[a], t[b], T_start_m, rising)
        sm = metrics(s["t"], s["Tm"], t[a], t[b], T_start_s, rising)
        fmt = lambda v: "never" if np.isnan(v) else f"{v:.0f}"
        print(f"  {name_j:10s} {tset[a]:5.1f}   {tset[a] - m[0]:+5.2f} / {tset[a] - sm[0]:+5.2f}"
              f"      {m[1]:5.1f} / {sm[1]:5.1f}          {m[2]:4.2f} / {sm[2]:4.2f}"
              f"            {m[3]:5.1f} / {sm[3]:5.1f}       {fmt(m[4]):>5s} / {fmt(sm[4]):>5s}")

    fig, (ax, axu, axe) = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.4, 1]})
    base = ~on & (t_all >= -BEFORE_S)
    ax.plot(t_all[base], run["temperature_C"][base], ".", color="0.75", ms=2.5,
            label="experiment: baseline (control off)")
    ax.plot(t, T, ".", color="0.6", ms=2.5, label=f"experiment: Run {n} (for comparison)")
    ax.plot(s["t"], s["Tm"], color="tab:blue", lw=1.6, label="simulation: measured mass $T_m$")
    ax.step(t, tset, "g--", where="post", lw=1.2, label="$T_{set}$")
    ax.axhline(t_amb, color="gray", ls=":", lw=1, label=f"$T_{{amb}}$ = {t_amb:.2f} °C (baseline)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"Run {n} simulation: {caption(c, kp, ki)}\n"
                 f"shared parameters: $G$ = {G_W_PER_C:g} W/K, $C_T/C$ = {CT_FRACTION:g}, "
                 f"$C$ = {c.total_capacitance_j_per_c:.1f} J/K, $H$ = {c.heat_loss_w_per_c:.2f} W/K, "
                 f"$P_{{u,h}}$ = {M.heating_tec_coefficient(c):.3f} W/PWM; $u_{{max}}$ = {u_max:g}, "
                 f"anti-windup {'on' if aw else 'off'}", fontsize=9.5)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.3)
    axu.step(t, u_sent, color="0.6", where="post", lw=1, label="experiment: PWM sent")
    axu.plot(s["t"], s["u_P"], color="tab:cyan", lw=1, label="simulation: $u_P$")
    if ki > 0:
        axu.plot(s["t"], s["u_I"], color="tab:orange", lw=1.4, label="simulation: $u_I$")
    axu.plot(s["t"], s["u"], color="k", lw=1.3, label="simulation: applied $u$ (signed)")
    axu.axhline(0, color="gray", lw=0.6)
    axu.set_ylabel("Signed PWM")
    axu.legend(fontsize=7.5, loc="upper right", ncol=2)
    axu.grid(alpha=0.3)
    axe.plot(t, run["error_C"][on], color="0.6", lw=1, label="experiment")
    axe.plot(s["t"], s["sp"] - s["Tm"], color="tab:blue", lw=1.2, label="simulation")
    axe.axhline(0, color="gray", lw=0.6)
    axe.set_ylabel("$e = T_{set} - T_m$ (°C)")
    axe.set_xlabel("Time since control started (s)")
    axe.legend(fontsize=7.5, loc="upper right")
    axe.grid(alpha=0.3)
    fig.tight_layout()
    out = MODULE6_FIG_DIR / f"part2_sim_run{n}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  Saved {out}")


def main():
    ap = argparse.ArgumentParser(description="Part II two-lump simulation Goldilocks plots")
    ap.add_argument("--run", type=int, choices=range(1, 7), help="one run (default: all six)")
    args = ap.parse_args()
    for n, name, _, _ in RUNS:
        if args.run in (None, n):
            do_run(n, name)


if __name__ == "__main__":
    main()
