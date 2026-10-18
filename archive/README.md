# Archive: programs not used in Part 7

Part 7 turns the Module 5 P controller into a PI controller on the hardware. The files
in this folder are not part of that program, so they were moved out of `python/` and
`arduino/` to keep the Part 7 work separate. Nothing was deleted and nothing was
changed except the paths: every script still runs from the repository root, and the
module note (`docs/module_notes/module_06_pi_control.md`) links here.

## Part 3 open-loop programs (superseded)

These were used once, for the open-loop steps of 17 Oct 2026
(`data/module_06/module4_run_20261017_164351.csv`, module note Part 3).

| File | What it is | Why it is archived |
| --- | --- | --- |
| `arduino/module4_tec_control/module4_tec_control.ino` | Module 4 sketch: `SET PWM <n> DIR HEAT/COOL` holds a fixed PWM; temperature shutdown, **no command timeout** | Open loop only. **Do not use it for P or PI control.** Without the timeout, the TEC keeps the last PWM if the GUI stops or the serial link drops. Part 7 uses `arduino/module6_pi_control/` (the Module 5 sketch), which has the timeout. |
| `python/module4_control_gui.py` | Module 4 manual-control GUI with the Module 5 `SerialLink` fix; logs to `data/module_06/openloop_<date>_<time>.csv` | Manual commands only, no feedback. The Part 7 GUI starts from the Module 5 P-control GUI instead; the serial link, reader thread, and line parser that Module 5 imported from this file are copied into `python/module6_pi_control_gui.py`, so nothing in `python/` imports it. |

To repeat an open-loop step: upload the sketch above, set `SERIAL_PORT` at the top of
the GUI, and run `python archive/python/module4_control_gui.py` (procedure in the note,
Part 3).

## Parts 3–6 analysis scripts (finished)

Each script makes the figures of one finished part of the note. They are kept so the
figures can be regenerated; rerunning them reproduces the committed figures.

| File | Note section | Input | Output in `docs/figures/module_06/` |
| --- | --- | --- | --- |
| `python/module6_tau63.py` | Part 3: 63 % estimate of $\tau$ with uncertainty | a Module 4 log or a Module 5 run (`--step`, `--t-final`, `--moving-start`, `--smooth`) | `tau63_<trace>[_step<n>].png` |
| `python/module6_part4_open_loop.py` | Part 4: measured open-loop steps vs the one-lump model | `data/module_06/module4_run_20261017_164351.csv` | `part4_open_loop_model_vs_measured.png` |
| `python/module6_part5_highgain.py` | Part 5: Module 5 $K_p$ = 32 trace vs the one-lump P model | `data/module_05/p_20261014-120316_Tset30.0_Kp32.csv` | `part5_highgain_Kp32_vs_one_lump.png` |
| `python/module6_part6_ki_sweep.py` | Part 6: $K_i$ sweep, $K_i = K_p/\tau$ pole–zero cancellation, Zero-integral recovery | none (one-lump PI model) | `part6_F_ki_sweep.png`, `part6_F_pole_zero.png` |
| `python/module6_part6_windup.py` | Part 6: windup with $u_{\max}$ = 30, anti-windup on and off | none (imports the model from `module6_part6_ki_sweep.py`) | `part6_G_windup.png` |

```bash
python archive/python/module6_tau63.py data/module_06/module4_run_20261017_164351.csv --step 1
python archive/python/module6_part4_open_loop.py    # imports module6_tau63.py
python archive/python/module6_part5_highgain.py
python archive/python/module6_part6_ki_sweep.py
python archive/python/module6_part6_windup.py       # imports module6_part6_ki_sweep.py
```

Two scripts import a neighbour, so keep the five scripts in the same folder. The full
list of `module6_tau63.py` commands used in the note is in Part 3.

## Not archived

- `python/Lab_6_7_modeling_tec.py`: the course simulation. It is still needed to
  choose $K_i$ before the Part 7 runs and for the Part II comparison, and the handout's
  commit command names this path.
- `data/`: the note links to every file there, and the Part 7 runs will be logged to
  `data/module_06/`.
- `docs/figures/module_06/`: the figures stay next to the note that shows them.
