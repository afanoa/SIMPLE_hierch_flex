---
name: rnasep-subclass-window-trajectory-2026-10-14
description: sliding-window trajectory over RNase P sub-classes along the synchronised 2D-latent coordinate (15 windows of 1/8 mass, stride 1/16, 42k particles each, 4-min cluster job): 8 disjoint windows reproduce at 0.31/0.38/0.20/0.15 per band, 2x the per-particle bins, about half the tracked states in the fine bands; the window family is 2-3 reproducible directions, not one; the mode-score ordering is unusable (data property, convention verified at 0.92)
metadata:
  type: project
---

Build: `mode_score.py` (window writer), `run_window_rec_sync.sbatch` (bin recipe with ARMS, K=15, odd windows copied
to a disjoint prefix for `group_maps_metric.py`), outputs `rnasep_coupled/trajectory/`, notes section 15 of
`doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md`, figures `viz_trajectory_sync_windows.png`,
`viz_trajectory_sync_corr.png`. Ordering = sub-class-level synchronised coordinate (near-only, 44-class component,
395 sub-classes, 329k particles; [[rnasep-sync-od-coordinate-2026-10-14]]).

Numbers (reproducible between-window variance, 117-30/30-20/20-12/12-8 A): windows 0.31/0.38/0.20/0.15; per-particle
bins of the same coordinate 0.18/0.16/0.08/0.10; tracked 15 states 0.36/0.58/0.38/0.23. Window family: first
eigenvalue 0.49 (coarse) to 0.29 (12-8 A); windows 1-2 and 10-11 carry opposite deviations of one arm, windows 14-15 a
second feature; window 11 has a localised strong deviation at 30-12 A (E/N 2.85 / 1.83), a candidate discrete
component worth a look. Mode-score ordering (sub-classes scored against the projected dominant mode of the tracked
family, `mode_volume.py` + `mode_score_abs.py`): even/odd reproducible 0.99 but tracked-state agreement r 0.15;
projection convention verified against the dumps' pm (0.92), so a single sub-class projection at 30-17 A is
dominated by other content than the one mode. Not used.

**Why:** pooling to sub-classes before windowing is what turns the noise-limited per-particle coordinate into a
reconstructable trajectory; the gain is the sub-class SNR.

**How to apply:** build trajectories at the sub-class level along a synchronised coordinate; use windows that
overlap for movies but disjoint windows for the metric (overlap inflates distance-1 correlation by construction);
the next lever is per-window pose refinement (poses determined per arm position), not a finer coordinate.
