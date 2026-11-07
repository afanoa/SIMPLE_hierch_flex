---
name: igm-arms-ceiling-pose-limited-2026-10-19
description: "IgM-like arms3: board E-step MATCHES the per-image ceiling (one round at true poses, 20 A: 43 % exact cell vs whitened-scorer 42 %); loss comes from refined poses (mean 17.6 deg, heavy tail), fine bands and EM contraction over rounds; arms are NOT below the SNR floor"
metadata:
  type: project
---

2026-10-19, arm_llr.py with fitted map scale + oracle_agree.py (board-numbering-correct readout):
- per-image ceiling at true poses (100 conformers): lp 10 d' 1.5, lp 15 1.8, lp 20 2.5; truth +2 deg d' 1.34, +5 deg 0.76.
- arms3 refined poses vs truth: direction median 7.5, mean 17.6 deg, 15 % > 20 deg, 8 % > 60 deg (55-conformer set:
  median 3.5, 1.3 % > 20). At refined poses lp 20 d' 0.96, lp 10 0.32.
- oracle region board, 7x3 lattice + 1 free, start at the true arm bins (exact-cell chance 0.048, within-1-step 0.30):
  true poses lp 20 ONE round: exact 0.43, within 1 step 0.85, corr(z1, theta bin) 0.86 = the whitened scorer on the
  same cells (0.42). Balance, free comp., pscale, |CTF| and focus each change < 0.03. Full 8-round schedule: exact 0.18
  (sigma 0.5-1.5) / 0.30 (sigma 0.5->0.2), corr 0.63 / 0.77: soft-EM contraction along the lattice. Refined poses,
  chain-like schedule to 8 A: exact 0.14, within 1 step 0.55, corr 0.54 (still informative, not random).
**Why:** earlier "drifted to random / below the SNR floor" claims came from an unscaled scorer and a label-numbering
bug ([[arm-llr-ceiling-needs-fitted-map-scale]], [[klem-oracle-labels-vs-board-numbering-trap]]).
**How to apply:** for hinged arms improve poses first (core-masked refinement), use a coarse band (~20 A) and narrow
final lattice sigmas; judge boards by corr(posterior mean, truth) or within-1-step agreement, not exact-cell purity.
