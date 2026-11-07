---
name: igm-motion-extremes-capture-2026-10-24
description: "IgM extremes + paths (2026-10-24): tree units miss arm extremes at >= 5k (0.41-0.49 of oracle); site-posterior tails 0.60-0.66; site-posterior quintile bins give monotone evenly stepped paths (best 0.70-0.75 of oracle span); tree-state paths coarse (latent) or wandering (map MST); latent tails weak"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-24T22:11:34.995Z
---

2026-10-24, phase 1 of the user's "are we capturing motion extremes" task, relayed by the coordinator; existing outputs only.
Tool: routes/toy55 `extremes.py` (+ extremes.sbatch); results in rnasep_klem/factorial/validation/extremes_igm.{txt,tsv} and
factorial_site_layer_results_2026_10_21.md.
- Metric: ratio = e(S) / e*(|S|), with |S| = S's genuine count. e(S) is the standardised mean of S's genuine members toward the end;
  e* is the mean of the |S| most extreme genuine particles.
- 20 ends (a1..a5 th / ph, low / high). The angles are near-uniform, so e* saturates at ~1.7 sd.
- Results (median best ratio):
  - leaves / units: m8c17 0.54 / 0.49 / 0.40 at n_min 2k / 5k / 10k. M: 0.59 at 2k, 0.41 at 5k, 0 at 10k (its units are mostly < 5k).
  - site posterior tails (top / bottom n by pole posterior): 0.60-0.66, with 3-7 ends >= 0.7.
  - site poles: 0.52 (M, 21 sites) / 0.62 (m8c17, 29 sites); the poles are 15k-120k halves.
  - root flex / fflex z tails: 0.35-0.46, with 26-49 % junk in the flex tails.
- a1_th low and a5_th high are weak everywhere; M's site layer misses a5_th high (0.08).
- PATHS (paths.py, same day; user: "how we go between them"):
  - Tree states as paths: the latent 15-NN graph gives 2-6 states with ~1.2-1.5 sd hops. The map MST gives 5.5-10 states, Spearman
    0.80-0.92, ~1 sd steps; M's wanders.
  - Continuous score: quintile bins of a kept site's posterior are monotone in most sites (18 of 21 M, 21 of 29 m8c17) and evenly
    stepped over ~2 sd; the best reach 0.70-0.75 of the oracle span.
  - Under empty-only, junk sits in the MIDDLE bins (up to 0.58).
  - The weak motions have deep candidates at abs(rho) 0.65-0.73 that the gate drops or that cover small levels: weak or narrow sites,
    not missing ones.
**Why:** the user wants the EXTREMES of each motion at good resolution ([[discrete-splits-capture-the-span-states-that-matter]]).
Tree partitions slice the middle of continuous motions, so the extremes need a continuous score per motion.
**How to apply:** for extremes, test per-site posterior top / bottom n (~10k) sets before tree units. Score any new method with
extremes.py against this baseline.
Related: [[m52-toy-arms-merged-units-2026-10-24]], [[toy-tree-validation-toolkit-2026-10-24]].
