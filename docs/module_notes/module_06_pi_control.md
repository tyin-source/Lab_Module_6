# Module 6, Part I: P/PI Control and Lumped Modeling

Handout: Phys 39 Fall 2026, Module 6, Part I (`labs/lab-06`). This note builds on
the Module 4 susceptibility and the Module 5 P-control runs. The data files used here
are copied unchanged from the [Module 5 repository](https://github.com/tyin-source/__Lab_Module_5__) into `data/module_04/` and
`data/module_05/`.

## Parameters carried forward

| Quantity | Value | Source |
| --- | --- | --- |
| $\chi_h$ (heating) | 0.485 °C/PWM count | Module 4 steady state vs PWM |
| $\chi_c$ (cooling) | 0.172 °C/PWM count | Module 4 steady state vs PWM |
| $r = \chi_h/\chi_c$ | 2.82 | Module 4 |
| $\tau$ (open loop) | ≈ 70 s (provisional, Part 3); to be measured from a dedicated open-loop step | Module 5 runs |
| $R_M$ (TEC module resistance) | 1.50 Ω at 27 °C | Laird data sheet (Module 4, A2) |
| $V_{\mathrm{TEC}}$ (on-state voltage across the TEC) | to be measured (Part 3 open-loop session) | |

**Sign convention.** As in the handout, $u$ is the signed PWM command: $u > 0$
heats and $u < 0$ cools. The open-loop relation $T = T_{\mathrm{amb}} + \chi u$ then
holds in both directions with the susceptibility of the active direction, so both
$\chi_h$ and $\chi_c$ are positive. For example, $u = -25$ gives
$T - T_{\mathrm{amb}} = 0.172 \times (-25) = -4.3$ °C. Module 5 sent an unsigned PWM
with a separate HEAT/COOL direction and therefore wrote the cooling value as
$\lvert\chi_c\rvert$. The number is the same.

Data: [`data/module_04/part4_steady_state.csv`](../../data/module_04/part4_steady_state.csv)
(Module 4 steady states) and [`data/module_05/p_runs_summary.csv`](../../data/module_05/p_runs_summary.csv)
(Module 5 run summary).

## Part 1: Algebraic droop model

### Derivation

Starting from the dimensional one-lump energy balance:

$$
C\frac{dT}{dt} = P_u u - H(T - T_{\mathrm{amb}})
$$

At steady state ($dT/dt = 0$), the open-loop temperature is:

$$
T = T_{\mathrm{amb}} + \chi u
$$

where the susceptibility is $\chi = P_u/H$, taken for the active direction
(this run heats, so $\chi = \chi_h$).

For P-only control, the command is $u = K_p(T_{\mathrm{set}} - T)$. Substituting
this into the open-loop relation:

$$
T = T_{\mathrm{amb}} + \chi K_p(T_{\mathrm{set}} - T)
$$

Collecting the terms in $T$ at steady state:

$$
(1 + \chi K_p)T_{\mathrm{ss}} = T_{\mathrm{amb}} + \chi K_p T_{\mathrm{set}}
$$

Subtracting from $T_{\mathrm{set}}$:

$$
\boxed{T_{\mathrm{set}} - T_{\mathrm{ss}} = \frac{T_{\mathrm{set}} - T_{\mathrm{amb}}}{1 + \chi K_p}}
$$

This is the droop formula, where the dimensionless loop gain is defined as
$L = \chi K_p$.

### Comparison with measurement

Run: Module 5, 14 Oct 2026, $K_p$ = 2 PWM/°C
([`p_20261014-113219_Tset30.0_Kp2.csv`](../../data/module_05/p_20261014-113219_Tset30.0_Kp2.csv)). $T_{\mathrm{amb}}$ is the
mean of the baseline before P control started; $T_{\mathrm{ss}}$ is the mean of
the last 120 s.

- $T_{\mathrm{set}}$ = 30.0 °C, $T_{\mathrm{amb}}$ = 22.94 °C,
  $\chi_h$ = 0.485 °C/PWM
- Loop gain $L = 0.485 \times 2 = 0.97$
- **Predicted** droop $= (30.0 - 22.94)/(1 + 0.97) = 3.58$ °C
- **Predicted** $T_{\mathrm{ss}}$ = 26.42 °C
- **Measured** droop = 3.76 °C ($T_{\mathrm{ss}}$ = 26.24 °C)

### Physical reason for the discrepancy

The measured droop is slightly larger than predicted (by 0.18 °C, 5 %). The
primary physical reason is that the baseline $T_{\mathrm{amb}}$ (22.94 °C) was
likely an overestimate. Prior to this run, the block was still relaxing toward a
cooled heat sink (baseline falling by 0.21 °C/min), meaning the effective
ambient temperature the block was seeing was lower than the baseline mean. A
lower $T_{\mathrm{amb}}$ would increase the numerator
$(T_{\mathrm{set}} - T_{\mathrm{amb}})$, thereby increasing the predicted droop
to match the measurement. Extrapolating the falling baseline with
$\tau$ ≈ 75 s gives $T_{\mathrm{amb}}$ ≈ 22.63 °C and a predicted droop of
3.74 °C, within 0.02 °C of the measured value (Module 5 note, Part 6 Q3).

## Part 2: Expressing the one-lump model using measured parameters

### Derivation of the measured parameters ($\chi$ and $\tau$)

We begin with the dimensional one-lump energy balance:

$$
C\frac{dT}{dt} = P_u u - H(T - T_{\mathrm{amb}})
$$

Dividing both sides by $C$, we get:

$$
\frac{dT}{dt} = \frac{P_u}{C}u - \frac{H}{C}(T - T_{\mathrm{amb}})
$$

At steady state, $dT/dt = 0$, which yields the steady-state temperature:

$$
T_{\mathrm{ss}} = T_{\mathrm{amb}} + \frac{P_u}{H}u
$$

Comparing this with the Module 4 open-loop relationship
$T = T_{\mathrm{amb}} + \chi u$, we identify the open-loop temperature
susceptibility:

$$
\boxed{\chi = \frac{P_u}{H}}
$$

To find the thermal time constant, we define the temperature displacement from
steady state as $\theta(t) = T(t) - T_{\mathrm{ss}}$. Since $T_{\mathrm{ss}}$ is
constant for a fixed command $u$, we have $d\theta/dt = dT/dt$. Substituting
$T = \theta + T_{\mathrm{ss}}$ into the energy balance, and using the
steady-state balance $P_u u = H(T_{\mathrm{ss}} - T_{\mathrm{amb}})$, we obtain:

$$
C\frac{d\theta}{dt} = -H\theta \quad\Longrightarrow\quad \frac{d\theta}{dt} = -\frac{H}{C}\theta
$$

The solution to this first-order differential equation is an exponential decay:

$$
\theta(t) = \theta(0)e^{-(H/C)t} = \theta(0)e^{-t/\tau}
$$

Comparing the exponents identifies the open-loop thermal time constant:

$$
\boxed{\tau = \frac{C}{H}}
$$

**Dimensional check:**

$$
[\tau] = \frac{\mathrm{J/K}}{\mathrm{W/K}} = \frac{\mathrm{J}}{\mathrm{J/s}} = \mathrm{s},
\qquad
[\chi] = \frac{\mathrm{W/PWM\ count}}{\mathrm{W/K}} = \frac{\mathrm{K}}{\mathrm{PWM\ count}} = \frac{^\circ\mathrm{C}}{\mathrm{PWM\ count}}
$$

### The measured-parameter model

The two measured quantities determine the two parameter ratios that appear in
the divided energy balance:

$$
\frac{H}{C} = \frac{1}{\tau},
\qquad
\frac{P_u}{C} = \frac{P_u/H}{C/H} = \frac{\chi}{\tau}
$$

Substituting them, the dimensional energy balance can be rewritten entirely in
terms of the measured quantities:

$$
\boxed{\frac{dT}{dt} = \frac{T_{\mathrm{amb}} + \chi u - T}{\tau}}
$$

This is not a different model; it is the one-lump energy balance expressed in a
form that can be simulated using only the measured susceptibility and time
constant. The quantity $(T_{\mathrm{amb}} + \chi u)$ is the steady temperature
toward which the model moves for a constant command $u$, and $\tau$ determines
how quickly it moves there.

With our measured values (Module 4 susceptibilities; $\tau$ from Part 3), using
the susceptibility for the active direction:

$$
\frac{dT}{dt} =
\begin{cases}
\dfrac{T_{\mathrm{amb}} + (0.485\ ^\circ\mathrm{C/PWM})\,u - T}{\tau}, & u > 0 \text{ (heating)}\\[2ex]
\dfrac{T_{\mathrm{amb}} + (0.172\ ^\circ\mathrm{C/PWM})\,u - T}{\tau}, & u < 0 \text{ (cooling)}
\end{cases}
\qquad \tau \approx 70\ \mathrm{s}
$$

Once $T_{\mathrm{amb}}$ and the command $u$ are specified, the model has no free
parameters, so its prediction can be compared directly with a measured trace
rather than fitted to it.

### Sidebar: dimensionless form

We can define dimensionless variables to generalize the model:

$$
\hat{t} = \frac{t}{\tau},
\qquad
\hat{T} = \frac{T - T_{\mathrm{amb}}}{\chi u_0},
\qquad
\hat{u} = \frac{u}{u_0}
$$

where $u_0$ is a reference command. Substituting these into the
measured-parameter model yields the dimensionless open-loop model:

$$
\frac{d\hat{T}}{d\hat{t}} = \hat{u} - \hat{T}
$$

This shows that $\tau$ supplies the time scale, and $\chi u_0$ supplies the
temperature scale. For P control, $u = K_p(T_{\mathrm{set}} - T)$ becomes
$\hat{u} = \chi K_p(\hat{T}_{\mathrm{set}} - \hat{T})$, so the only remaining
parameter is the dimensionless loop gain:

$$
L = \chi K_p
$$

## Part 3: Estimate $\tau$

**Method (63 %).** For a step in the applied PWM,
$T_{63} = T_{\mathrm{initial}} + 0.63\,(T_{\mathrm{final}} - T_{\mathrm{initial}})$, and
$\tau$ is the time from the step until the temperature first reaches $T_{63}$.
`python/module6_tau63.py` does this for a trace. $T_{\mathrm{initial}}$ is the mean of
the 30 s before the step and $T_{\mathrm{final}}$ is the mean of the last 60 s. The
crossing is taken on a 5-reading centred running mean, so that a single noisy reading
does not cross early. The uncertainty is the noise and the $T_{\mathrm{final}}$
uncertainty divided by the slope at $T_{63}$. Every reading in the Module 5 files is
an average of 1,000 thermistor-voltage samples.

**Available traces.** Module 4 saved only the steady temperatures, not the time
series, and the Module 5 GUI logged only P-control runs. No saved trace is a large
open-loop step that was held until it settled. Three estimates are possible from the
Module 5 files:

1. **Closed-loop P steps (main estimate).** In the one-lump model a P step relaxes
   with $\tau_{\mathrm{cl}} = \tau/(1+L)$ (Part 5), so
   $\tau = t_{63}(1+L)$ with $L = \chi_h K_p$. Only runs that started from a settled
   block and did not overshoot are used:

   | $K_p$ (PWM/°C) | $L$ | $t_{63}$ (s) | $\tau = t_{63}(1+L)$ (s) | Run |
   | ---: | ---: | ---: | ---: | --- |
   | 1 | 0.49 | 48.4 | 72 | [`p_20261017-005209_Tset30.0_Kp1.csv`](../../data/module_05/p_20261017-005209_Tset30.0_Kp1.csv) |
   | 2 | 0.97 | 35.5 | 70 | [`p_20261014-113219_Tset30.0_Kp2.csv`](../../data/module_05/p_20261014-113219_Tset30.0_Kp2.csv) |
   | 4 | 1.94 | 22.6 | 66 | [`p_20261014-114135_Tset30.0_Kp4.csv`](../../data/module_05/p_20261014-114135_Tset30.0_Kp4.csv) |
   | 8 | 3.88 | 15.1 | 74 | [`p_20261014-114820_Tset30.0_Kp8.csv`](../../data/module_05/p_20261014-114820_Tset30.0_Kp8.csv) |

   Mean $\tau$ = 70 s, spread ±4 s. $K_p$ = 0.5 is left out because its PWM fell
   from 4 to 3 after 54 s, so it is not a single step. $K_p$ = 16 and 32 are left out
   because they overshoot.
2. **Open-loop cooling step, PWM 0 → −2** (sign test,
   [`p_20261014-110522_Tset18.0_Kp0.5.csv`](../../data/module_05/p_20261014-110522_Tset18.0_Kp0.5.csv)). The applied PWM stayed at −2, so this is
   a true open-loop step, but it is only 0.34 °C deep and lasted 131 s. With
   $T_{\mathrm{final}} = T_{\mathrm{amb}} - \chi_c \cdot 2 = 22.33$ °C, the 63 % time
   is $\tau \approx 80 \pm 20$ s. The trace is still falling almost linearly at
   131 s and is already below 22.33 °C, so it had not settled and $T_{\mathrm{final}}$
   is uncertain.
3. **Constant PWM 3 after 54 s** in the 17 Oct $K_p$ = 0.5 rerun,
   [`p_20261017-003102_Tset30.0_Kp0.5.csv`](../../data/module_05/p_20261017-003102_Tset30.0_Kp0.5.csv) (open loop in
   practice, because rounding held the command at 3; $T_{\mathrm{initial}}$ is taken
   at the switch from PWM 4 to 3). The remaining approach is only
   0.45 °C, and the room wander of about ±0.1 °C at the end gives
   $\tau \approx 87 \pm 56$ s.

![63 % method on the PWM −2 cooling step](../figures/module_06/tau63_p_20261014-110522_Tset18.0_Kp0.5.png)

**Provisional result:** $\tau \approx 70$ s. The uncertainty is about ±10 %:

- **Spread:** ±4 s between the four runs.
- **Loop delay:** the delay of about 3 readings at the start of each run adds about
  2 s to $t_{63}$, which becomes $2(1+L)$ s in $\tau$. That is up to about 10 s at
  $K_p$ = 8, so these values are upper estimates.
- **Low-PWM slope:** the slope below PWM 12 is about 7 % smaller than the fitted
  0.485, which changes $1+L$ by 3–6 %.

The estimate is model-dependent because it assumes $\tau_{\mathrm{cl}} = \tau/(1+L)$.
The fitted closed-loop time constants of the same runs give
$\tau_{\mathrm{cl}}(1+L)$ = 69–77 s (Module 5, Part 6). As the handout notes, this
fixes only the ratio $C/H$, not $C$ and $H$ separately.

Reproduce (from the repository root):

```bash
python python/module6_tau63.py data/module_05/p_20261014-110522_Tset18.0_Kp0.5.csv --t-final 22.333 --t-final-err 0.05
python python/module6_tau63.py data/module_05/p_20261017-003102_Tset30.0_Kp0.5.csv --step 2 --moving-start
```

The closed-loop $t_{63}$ values are the `tau63_s` column of
[`data/module_05/p_runs_summary.csv`](../../data/module_05/p_runs_summary.csv).

### Open-loop step measurement (to do)

A large open-loop step gives $\tau$ directly, without $\chi$, $K_p$, or
$\tau_{\mathrm{cl}} = \tau/(1+L)$. The same session supplies the measured open-loop
heating and cooling traces that Part 4 compares with the simulation, and the on-state
TEC voltage.

**Tools.** Use the Module 4 sketch `arduino/module4_tec_control/` with
`python/module4_control_gui.py`. The Module 5 sketch cannot be used here: it stops the
TEC when no command arrives for 3 s, and this GUI sends a command only when the PWM or
direction changes. The GUI logs every reading (a 1,000-sample average) to
`data/module_06/openloop_<date>_<time>.csv`.

**Before the steps.**

1. Upload `arduino/module4_tec_control/module4_tec_control.ino`. The start-up banner
   must show the 60 °C software limit.
2. Set `SERIAL_PORT` at the top of `python/module4_control_gui.py` and close the
   Arduino Serial Monitor.
3. Leave the block at PWM 0 for at least 15 min after any heating run.
4. From the repository root, run `python python/module4_control_gui.py`. Logging
   starts at once. Record at least 2 min at PWM 0; the temperature should drift by
   less than about 0.05 °C/min. Note the room temperature and the clock time.

**Steps.** Type the PWM in the text box and press Enter, so that each change is a
single command. Change the direction only while the PWM is 0.

| Step | Command | Hold | Expected |
| ---: | --- | --- | --- |
| 1 | HEAT, PWM 0 → 25 | ≥ 8 min (about 7 τ) | rise $\chi_h \cdot 25 \approx 12$ °C, to about 34 °C |
| 2 | PWM 25 → 0 | ≥ 15 min | back to room temperature; a heat-sink dip below room temperature is possible |
| 3 | select COOL at PWM 0, then PWM 0 → 25 | ≥ 8 min | fall $\chi_c \cdot 25 \approx 4.3$ °C, to about 17.7 °C |
| 4 | PWM 25 → 0 | ≥ 10 min | back to room temperature |

5. **On-state TEC voltage.** With the block back near room temperature, set COOL,
   PWM 255 for at most 5 s, read the DC voltage across the TEC terminals with the
   multimeter, then set PWM 0. At a lower PWM the multimeter reads the duty-cycle
   average, not the on-state voltage. Expect the block to cool by about 3 °C during
   the 5 s.
6. Close the GUI and record the file name, room temperature, and $V_{\mathrm{TEC}}$.
   Upload the Module 5 sketch again before any P or PI run, because it has the
   command timeout.

Set PWM 0 at once if the temperature moves the wrong way, goes above 45 °C, the
display freezes, or a safety shutdown is reported.

**Analysis.**

```bash
python python/module6_tau63.py data/module_06/openloop_<date>_<time>.csv --step 1   # 0 -> +25
python python/module6_tau63.py data/module_06/openloop_<date>_<time>.csv --step 3   # 0 -> -25
python python/module6_tau63.py data/module_06/openloop_<date>_<time>.csv --step 2   # +25 -> 0
python python/module6_tau63.py data/module_06/openloop_<date>_<time>.csv --step 4   # -25 -> 0
```

Steps 1 and 3 start from a settled block and are the main results. Steps 2 and 4 are
relaxations at PWM 0; the heat-sink dip can bias their $T_{\mathrm{final}}$. The
voltage pulse adds two more steps at the end, which are not used. For a 12 °C step,
noise moves the 63 % time by well under 1 s. Whether heating and cooling give the
same $\tau$ is itself a test of the one-lump model, in which $\tau = C/H$ does not
depend on direction.
