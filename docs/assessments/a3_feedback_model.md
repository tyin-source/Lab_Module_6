# A3: Feedback Data and Lumped-Model Memo

Tianren Yin · Phys 39, Fall 2026 · due Wednesday 21 October, 6:00 PM (Moodle PDF
`A3_Yin.pdf`).

The memo is written in LaTeX (REVTeX 4.2, as A2):

| File | What it is |
| --- | --- |
| [`a3/A3.tex`](a3/A3.tex) | Source of the memo: Questions 1–5, data and code links |
| [`a3/A3.pdf`](a3/A3.pdf) | Compiled memo (rename to `A3_Yin.pdf` for Moodle) |
| [`a3/make_figures.py`](a3/make_figures.py) | Makes the three figures in `a3/figures/` |

Build (from the repository root):

```bash
python docs/assessments/a3/make_figures.py      # figures
cd docs/assessments/a3 && latexmk -pdf A3.tex   # PDF (Overleaf also works: upload a3/)
```

## Where each answer's evidence lives

| Question | Evidence | Code |
| --- | --- | --- |
| 1. One-lump model, P droop, $\tau_{\mathrm{cl}}$ | Part I note, Parts 1–3; Table I | `python/module6_a3_metrics.py` |
| 2. Six experimental Goldilocks plots, step metrics | Part I note, Part 7; Fig. 1, Table II | `python/module6_part7_runs.py`, `python/module6_a3_metrics.py` |
| 3. Six two-lump simulation plots, comparison | Part II note, Part 2; Fig. 2, Table III, Fig. 3 | `python/module6_part2_goldilocks.py` |
| 4. PI steady state, windup, anti-windup | Part I note, Part 6 (integral, windup) and Part 7 (program) | `python/module6_pi_control_gui.py` (`pi_control`) |
| 5. Model limits, gain choice, loop delay | Part II note, Steps 2–5; Fig. 3 | `python/module6_a3_delay_check.py` |

Notes: [`docs/module_notes/module_06_pi_control.md`](../module_notes/module_06_pi_control.md)
(Part I) and [`docs/module_notes/module_06_part2_two_lump.md`](../module_notes/module_06_part2_two_lump.md)
(Part II). Run files: [`data/module_06/`](../../data/module_06/) (the six retained runs are
listed in its README). Modeling checkpoint cited in the memo: commit `0ebd3a5`.
