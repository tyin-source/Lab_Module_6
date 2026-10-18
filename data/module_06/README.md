# Module 6 measurements

| File | What it is | Used as |
| --- | --- | --- |
| `module4_run_20261017_164351.csv` | Part 3 open-loop steps (0 → +25 → 0 → −25 → 0 PWM), 17 Oct 2026 | Part 3 $\tau$, Part 4 model check |
| `p_20261018-045404_Kp2_Ki0_T30-35.csv` | Part 7 P, $K_p$ = 2: 30 → 35 → 30 °C | **Goldilocks Run 1** |
| `p_20261018-054301_Kp8_Ki0_T30-35.csv` | Part 7 P, $K_p$ = 8: 30 → 35 → 30 °C | **Goldilocks Run 2** |
| `p_20261018-062146_Kp32_Ki0_T30-35.csv` | Part 7 P, $K_p$ = 32: 30 → 35 → 30 °C | **Goldilocks Run 3** |
| `pi_20261018-063831_Kp6_Ki0.04_T30-35.csv` | Part 7 PI, $K_p$ = 6, $K_i$ = 0.04: 30 → 35 → 30 °C | **Goldilocks Run 4** |
| `pi_20261018-072600_Kp6_Ki0.08_T30-35.csv` | Part 7 PI, $K_p$ = 6, $K_i$ = 0.08: 30 → 35 → 30 °C | **Goldilocks Run 5** |
| `pi_20261018-083656_Kp6_Ki0.5_T30-35.csv` | Part 7 PI, $K_p$ = 6, $K_i$ = 0.5: 30 → 35 → 30 °C | **Goldilocks Run 6** |
| `p_20261018-034128_Kp2_Ki0_T25-30.csv` | $K_p$ = 2, 25 °C setpoint, stopped after 30 s | debugging only |
| `p_20261018-034247_Kp2_Ki0_T25-30.csv` | $K_p$ = 2, 25 °C setpoint only; setpoints then changed to 30/35 °C | debugging only |
| `p_20261018-035802_Kp2_Ki0_T30-35.csv` | $K_p$ = 2, 30 °C setpoint only | debugging only |
| `p_20261018-041612_Kp2_Ki0_T30-35.csv` | $K_p$ = 2, started at 35 °C from room temperature (no 30 → 35 step) | debugging only |
| `p_20261018-043355_Kp2_Ki0_T30-35.csv` | $K_p$ = 2, 30 → 35 °C, stopped before the downward step | debugging only |
| `p_20261018-061226_Kp32_Ki0_T30-35.csv` | $K_p$ = 32, stopped after 1 s | debugging only |
| `pi_20261018-075901_Kp6_Ki0.05_T30-35.csv` | Intended as Run 6 ($K_i$ = 0.50) but run with $K_i$ = 0.05: slow, no overshoot, between Runs 4 and 5 | debugging only |
| `p_20261018-051855_Kp8_Ki0_T28-32.csv` | $K_p$ = 8, 28 → 32 → 28 °C (setpoints differ from the study); repeated at 30/35 °C | debugging only |

Run files from `python/module6_pi_control_gui.py` are named
`<p|pi>_<date>-<time>_Kp<..>_Ki<..>_T<low>-<high>.csv`; the columns are listed at the top
of that file. Debugging runs did not meet the matched conditions of the Goldilocks
study and are not used in A3 (handout, Part 7).
