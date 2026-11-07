---
name: flex-pcg-falloff-prior-pulls-solvent-shell-2026-10-22
description: "pcg_falloff alone on RNase P r1 pulls a SOLVENT shell, wrecks latent/poles; + pcg_complement: reproducible region latent, but held out its closed pole is LESS closed than hard's (imprint); keep the hard support"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T15:04:50.037Z
---

User 2026-10-22: "the flex only has a solvent support ... it needs a falloff prior like in the refine3D". I implemented it as
`pcg_falloff`:
- branch flex-falloff-prior, uncommitted;
- write-up in doc/implementation_notes/flex_pca_falloff_prior.md §6-7.

I tested it on node r1's fflex (batch2 stree_port/falloff_r1): hard support vs prior, same inputs, same node mask.

Results:
- **The halves do agree on outside-mask variance.** CV picks lambda_rel 0.03-0.13, and a hard cut costs 68 % in cross-half error.
- **But that variance is a solvent shell.** Of the eigenvalue-weighted basis energy:
  - 0.647 is below the mean's Otsu level, peaking at 75-90 A from the centre;
  - 0.018 is molecule outside the mask;
  - 0.000 is on the closed-arm region or the tip.
- **The shell takes over the latent.** Latent CCA hard vs prior runs from 0.64 down to 0.03. Heterogeneity resolution drops from 7.3 to
  ~20 A. The poles stop reproducing p2: the closed pole goes from 13,377 particles to a 947-particle core, 98 % of which lies inside p2's
  closed pole.
- **The closed pole's tip coverage at 95-110 A is unchanged:** 0.71/0.67/0.59 hard vs 0.74/0.71/0.58 prior. At the fflex box the r1
  node mask already covers the closed arm (mean w 0.85, tip 0.38, nothing below the 0.05 floor).
- **The noise scale stays put:** 0.41 -> 0.44. Do NOT read the fits' sig2 as a residual: sig2_eff is the images' high-frequency noise
  floor (em_calibrate_noise_prior, power above 0.7 Nyquist), never refitted from the model, so no basis change can move it
  [[flex-pca-deconv-precision-mismatch-flag-2026-10-22]].

**Why:** cross-half CV separates signal from noise, not relevant signal from irrelevant. Reproducible solvent-shell variance passes it
and swamps the region's components.

**Complement result (same day, `pcg_complement=yes`, ported on flex-falloff-prior):**
- The complement block takes 66 % of the basis energy, 3/4 solvent; the region block goes back inside the mask (93 % core + tail).
- The region latent is as reproducible across seeds as hard (complement seed 1 vs 2: CCA 0.997-0.990; hard vs port2 0.998-0.988),
  but differs systematically from hard (CCA 0.84 -> 0.38 on BOTH seeds). Read it as outside variance removed from the region latent,
  not as noise.
- Poles reproduce across seeds (J 0.88 / 0.81 / 0.71); the closed pole shrinks to 2.5-3.4k (J 0.62) and reads more closed (tip
  0.67-0.73 at 105-110 A vs 0.59), but those are selected-subset maps [[selected-subset-maps-imprint-inflates-extremes-2026-10-22]].
- Noise scale 0.37; deconv K unstable (5 vs 2).

**How to apply:**
- Never run pcg_falloff without pcg_complement.
- Judge the complement against hard by CROSS-SEED reproducibility and held-out extremes, never by CCA against hard: a cleaner latent
  also lowers that CCA.
- HELD-OUT VERDICT (pixel-split, lsamp2h_poleheld.py): the complement's closed pole is LESS closed than hard's at matched N
  (0.66-0.70 / open 0.67-0.73 vs 0.74-0.75 / 0.54-0.67 at N 3,000); only 15-21 % of a held-out re-selection falls in its pole.
  Both reproduce across seeds within noise (cross-fold rho 0.88-1.00); hard vs complement differ beyond noise. KEEP THE HARD SUPPORT.
- Compare held-out maps ACROSS pixel-set folds (B-built vs A-built); within-fold comparisons share noise and give rho > 1.
- Before blaming the support, check the mask's coverage of the moving region at the fflex box.
- STATUS (user, 2026-10-22): both switches kept as off-by-default options, uncommitted on flex-falloff-prior. Open item: a
  membrane-protein test (cnga1 needs a re-import first [[pre-release4-projects-fail-flex-contracts-2026-10-22]]; scripts ready
  in cnga1/falloff_cmpl).

Related: [[focus-mask-truncates-mobile-tip-falloff-prior-2026-10-21]], the mask-building falloff for pcg_focus. That is a different
mechanism, and it still stands.
