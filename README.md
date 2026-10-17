# Phys 39 Lab: TEC Temperature Control (Module 6)

Module 6, Part I develops the one-lump model of the TEC block and adds PI control;
Part II compares one- and two-lump models with the course simulation.

## Contents

| Path | What it is |
| --- | --- |
| `docs/module_notes/module_06_pi_control.md` | Module note (Part I, Parts 1–3 so far) |
| `python/module6_tau63.py` | 63 % time-constant estimate from a PWM step, with uncertainty |
| `python/module4_control_gui.py` | Module 4 open-loop GUI (with the Module 5 `SerialLink` fix); logs to `data/module_06/openloop_<date>_<time>.csv` |
| `arduino/module4_tec_control/` | Module 4 sketch (no command timeout) for open-loop steps |
| `data/module_04/`, `data/module_05/` | Module 4 and Module 5 data used here, copied from the Module 5 repository |
| `data/module_06/` | Module 6 measurements |
| `docs/figures/module_06/` | Figures |

## How to run

```bash
pip install -r requirements.txt
python python/module4_control_gui.py          # open-loop steps (Module 4 sketch)
python python/module6_tau63.py <trace.csv> --step 1
```

Set `SERIAL_PORT` at the top of the GUI first. The open-loop step procedure is in the
module note, Part 3.
