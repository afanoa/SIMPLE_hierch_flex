---
name: igm-focused-flex-mask-is-bottleneck-2026-11-07
description: "VERDICT T1 2026-11-07: on the IgM fixture a flex_pca fit focused on the TRUE arm mask sees all five arms (R2 .19-.31, LDA NMI .11-.21) where the whole-molecule latent reads 0; the tree's island masks (13-17 % on the arm, rest solvent) fail -> site discovery/masks are the bottleneck"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T07:12:42.135Z
---

T1 of guide section 7.7 (fork hetero-tree-wip), run 2026-11-07 on `igm_cryolo/tree/fixture_igm` (5 x 50 kDa arms, discrete hinges,
null twin continuous). Files: `fixture_igm/t1_sitefocus/T1_VERDICT.md`; generic scripts `t1_sitefocus/scripts/sitefocus_run.sbatch`,
`sitefocus_score.py` (`--nodecoords` adds node + site), `t1_prep.py`, `t1_table.py`.
- True arm mask (`pcg_mskfile`), root poses, the tree's own flex keys: every arm read, R2 .19-.31, cross-validated LDA NMI a1 .211,
  a2-a5 .11-.14 (cross-validated site stage B: a1 .146, a2 .047, a5 .039, a3/a4 none); unsupervised best-axis split a1 .144, a3 .119,
  a5 .108, a2 .104, a4 .006; no leakage between arms. Whole-molecule latent: 0.000 on every arm.
- The tree's island masks: R2 .04-.11, LDA NMI .016-.063; each island 13-17 % on its arm, rest solvent; heavy-tailed coordinates
  (kurtosis 5-28) in both twins.
- Per-particle discreteness unreadable at this SNR (tied-variance 2-vs-1 gain ~0 for discrete and continuous arms); free-variance
  gains are bought by heavy tails (same trap as the sitefocus-gates session (`hetero_tree_20261028/handoff/sitefocus-gates/VERDICT.md`) found on qtoynull). Discreteness stays map-level.
- Gate fixtures (sitefocus-gates): focus does not beat the node latent on r131/qtoy2; node + sites never worse; qtoyflex gains
  (closed .786 -> .845).

**Why:** confirms the user's thesis that site/mask is the weakest layer: the latent can read small moving parts only through
arm-shaped masks.
**How to apply:** next = data-derived arm-shaped masks (islands trimmed to the molecule, or refined from the focused fit's mode
energy: T1b, `fixture_igm/t1b_masks/`), then focused coordinates appended to the node's. Never use free-variance gains to call
discreteness or rank sites. Related: [[site-mask-rework-direction-2026-11-07]], [[igm-fixture-xval-stageb-gapw-2026-11-06]].
