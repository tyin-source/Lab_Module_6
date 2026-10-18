# Phys 39 Module 6, Part 3: estimate the thermal time constant with the 63 % method.
#
# For a step in the applied PWM:
#   T_initial = mean of the readings in the 30 s before the step (after any earlier
#               change); with --moving-start, the readings just around the step, for a
#               step taken before the previous level had settled
#   T_final   = mean of the last FINAL_WINDOW_S of the trace (or --t-final if the
#               trace was stopped before it settled)
#   T_63      = T_initial + 0.63 (T_final - T_initial)
#   tau       = time from the step to the first reading at or past T_63
#               (linearly interpolated between readings)
#
# For a small step the noise makes single readings cross T_63 early, so the crossing is
# taken on a centred running mean of --smooth readings (default 5); the raw first
# crossing is printed as well. A centred mean does not shift a smooth trace in time.
#
# Uncertainty: the noise of the (smoothed) readings and the uncertainty in T_final both
# move the crossing time by (temperature error) / (slope of the trace at T_63).
#
# Reads either a Module 4 GUI log (time_s, temperature_C, pwm, heat_cool, ...) or a
# Module 5 run file (time_s, temperature_C, ..., pwm_arduino, heat_cool, ...). The
# step is the n-th change of the applied signed PWM (--step, default the first).
#
# Usage (from the repository root):
#   python python/module6_tau63.py <trace.csv> [--step 2] [--t-final 22.33 --t-final-err 0.05]
#   Open-loop session (data/module_06/openloop_*.csv): --step 1 is 0 -> +25, --step 2 is
#   +25 -> 0, --step 3 is 0 -> -25, --step 4 is -25 -> 0.
# Figure: docs/figures/module_06/tau63_<trace>.png

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = REPO_ROOT / "docs" / "figures" / "module_06"
BASELINE_S = 30.0      # readings before the step used for T_initial and the noise
FINAL_WINDOW_S = 60.0  # end of the trace used for T_final


def read_trace(path: Path):
    """Return time (s), temperature (C) and applied signed PWM from either CSV format."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r["time_s"]) for r in rows])
    T = np.array([float(r["temperature_C"]) for r in rows])
    pwm_key = "pwm_arduino" if "pwm_arduino" in rows[0] else "pwm"
    pwm = np.array([float(r[pwm_key]) for r in rows])
    heat = np.array([int(float(r["heat_cool"])) == 1 for r in rows])
    return t, T, np.where(heat, pwm, -pwm)


def first_crossing(t, T, level, rising):
    """Interpolated time of the first reading at or past level."""
    past = T >= level if rising else T <= level
    k = int(np.argmax(past))
    if not past[k]:
        return np.nan
    if k == 0:
        return t[0]
    f = (level - T[k - 1]) / (T[k] - T[k - 1])
    return t[k - 1] + f * (t[k] - t[k - 1])


def main():
    ap = argparse.ArgumentParser(description="63 % time-constant estimate from a PWM step")
    ap.add_argument("trace", type=Path)
    ap.add_argument("--step", type=int, default=1, help="use the n-th change of applied PWM")
    ap.add_argument("--t-final", type=float, default=None,
                    help="final temperature (C) if the trace did not settle")
    ap.add_argument("--t-final-err", type=float, default=None,
                    help="uncertainty of T_final (C); default: noise and drift of the final window")
    ap.add_argument("--moving-start", action="store_true",
                    help="previous level had not settled: T_initial from the readings at the step")
    ap.add_argument("--smooth", type=int, default=5,
                    help="centred running mean (readings) used for the crossing; 1 = raw")
    args = ap.parse_args()

    t, T, u = read_trace(args.trace)
    changes = np.where(np.diff(u) != 0)[0] + 1
    if len(changes) < args.step:
        raise SystemExit(f"only {len(changes)} PWM change(s) in {args.trace.name}")
    k0 = changes[args.step - 1]
    t_step = t[k0]
    nxt = changes[args.step] if len(changes) > args.step else len(t)
    u_before, u_after = u[k0 - 1], u[k0]
    held_s = t[nxt - 1] - t_step

    t_prev = t[changes[args.step - 2]] if args.step > 1 else -np.inf
    before = (t < t_step) & (t >= max(t_step - BASELINE_S, t_prev))
    base = (t < t[changes[0]]) & (t >= t[changes[0]] - BASELINE_S)
    tb = t[base]
    noise = float(np.std(T[base] - np.polyval(np.polyfit(tb, T[base], 1), tb)))
    if args.moving_start or before.sum() < 5:
        T_i = float(T[k0 - 2:k0 + 1].mean())          # previous level was still moving
    else:
        T_i = float(T[before].mean())

    seg_t, seg_T = t[k0:nxt] - t_step, T[k0:nxt]
    final = seg_t >= seg_t[-1] - FINAL_WINDOW_S
    if args.t_final is not None:
        T_f = args.t_final
        T_f_err = args.t_final_err if args.t_final_err is not None else 0.0
    else:
        T_f = float(seg_T[final].mean())
        drift = abs(np.polyfit(seg_t[final], seg_T[final], 1)[0]) * FINAL_WINDOW_S
        T_f_err = args.t_final_err if args.t_final_err is not None else max(
            float(seg_T[final].std()), drift)

    rising = T_f > T_i
    T63 = T_i + 0.63 * (T_f - T_i)
    n = max(1, args.smooth)
    smooth = np.convolve(seg_T, np.ones(n) / n, mode="same")
    smooth[: n // 2] = seg_T[: n // 2]                # no partial windows at the ends
    if n // 2:
        smooth[-(n // 2):] = seg_T[-(n // 2):]
    tau_raw = first_crossing(seg_t, seg_T, T63, rising)
    tau = first_crossing(seg_t, smooth, T63, rising)
    slope = abs(T_f - T_i) * 0.37 / tau if np.isfinite(tau) else np.nan
    err_noise = noise / np.sqrt(n) / slope
    err_final = 0.63 * T_f_err / slope
    err = float(np.hypot(err_noise, err_final))

    print(f"{args.trace.name}: step {args.step}, PWM {u_before:+.0f} -> {u_after:+.0f} at "
          f"t = {t_step:.1f} s, held {held_s:.0f} s")
    print(f"  T_initial {T_i:.3f} C (noise {noise:.3f} C), T_final {T_f:.3f} +/- {T_f_err:.3f} C, "
          f"change {T_f - T_i:+.3f} C")
    print(f"  T_63 {T63:.3f} C first reached {tau:.1f} s after the step "
          f"({n}-reading mean; raw readings {tau_raw:.1f} s)")
    print(f"  tau = {tau:.1f} +/- {err:.1f} s  (noise {err_noise:.1f} s, T_final {err_final:.1f} s)")
    if held_s < 3 * tau:
        print(f"  note: PWM was held only {held_s / tau:.1f} tau; T_final is not a settled value")

    fig, (ax, axp) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    rel = t - t_step
    ax.plot(rel, T, "k.", ms=3, label="measured")
    if n > 1:
        ax.plot(seg_t, smooth, "tab:blue", lw=1, label=f"{n}-reading running mean")
    ax.axhline(T_i, color="gray", ls=":", lw=1, label=f"T_initial = {T_i:.2f} °C")
    ax.axhline(T_f, color="gray", ls="--", lw=1, label=f"T_final = {T_f:.2f} °C")
    ax.axhline(T63, color="tab:red", ls="-", lw=1, label=f"T_63 = {T63:.2f} °C")
    ax.axvline(tau, color="tab:red", ls="-", lw=1)
    err_txt = f"{err:.1f}" if err < 1 else f"{err:.0f}"
    ax.annotate(f"τ ≈ {tau:.0f} ± {err_txt} s", (tau, T63), xytext=(8, -14),
                textcoords="offset points", color="tab:red")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"63 % method, step {u_before:+.0f} → {u_after:+.0f} PWM\n{args.trace.name}",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    axp.step(rel, u, "k", where="post", lw=1)
    axp.set_xlim(-BASELINE_S, held_s + 10)            # this step only
    shown = (rel >= -BASELINE_S) & (rel <= held_s + 10)
    pad = 0.05 * (T[shown].max() - T[shown].min()) + 0.05
    ax.set_ylim(T[shown].min() - pad, T[shown].max() + pad)
    axp.set_ylabel("Applied signed PWM")
    axp.set_xlabel("Time since the step (s)")
    axp.grid(alpha=0.3)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / f"tau63_{args.trace.stem}{'' if args.step == 1 else f'_step{args.step}'}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
