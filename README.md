# Phys 39 Lab: TEC Temperature Control (Module 6)

Module 6, Part I develops the one-lump model of the TEC block and adds PI control;
Part II compares one- and two-lump models with the course simulation.

## Contents

| Path | What it is |
| --- | --- |
| `docs/module_notes/module_06_pi_control.md` | Module note, Part I (Parts 1–7) |
| `docs/module_notes/module_06_part2_two_lump.md` | Module note, Part II (two-lump model and simulation) |
| `python/module6_part2_lienhard_1_8.py` | Part II pre-class: Lienhard Problem 1.8 (lumped cooling of a copper sphere) |
| `python/module6_pi_control_gui.py` | Part 7 P/PI control GUI (from the Module 5 P-control GUI); logs to `data/module_06/` |
| `python/module6_common.py` | Module 4 susceptibilities, $\tau$, one-lump PI damping ratio, run-file reader |
| `python/module6_pi_gui_selftest.py` | Tests the GUI without hardware against a simulated one-lump Arduino |
| `arduino/module6_pi_control/` | Module 5 sketch for P/PI runs (temperature shutdown, 0–255 clamp, 3 s command timeout) |
| `python/Lab_6_7_modeling_tec.py` | Course simulation (one- and two-lump models, open-loop, P, and PI control) |
| `archive/` | Programs not used in Part 7: the Part 3 open-loop sketch and GUI, and the Parts 3–6 analysis scripts (see `archive/README.md`) |
| `data/module_04/`, `data/module_05/` | Module 4 and Module 5 data used here, copied from the Module 5 repository |
| `data/module_06/` | Module 6 measurements |
| `docs/figures/module_06/` | Figures |


## How to run

```bash
pip install -r requirements.txt
python python/module6_pi_gui_selftest.py             # check the P/PI GUI (no hardware)
python python/module6_pi_control_gui.py              # P/PI runs (Module 6 sketch uploaded)
python python/Lab_6_7_modeling_tec.py                # course simulation
python archive/python/module6_part6_ki_sweep.py      # e.g. regenerate a Part 6 figure
```

Set `SERIAL_PORT` at the top of the GUI first. The Part 7 procedure is in the module
note, Part 7. The commands for every archived script, and the open-loop step
procedure, are in `archive/README.md` and the module note, Part 3.
