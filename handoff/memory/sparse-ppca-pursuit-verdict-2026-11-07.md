---
name: sparse-ppca-pursuit-verdict-2026-11-07
description: "VERDICT YES (lever 6, guide 7.8): group-sparse PPCA by greedy site pursuit learns IgM arm masks + coordinates from data (10 % fix: 5/5 arms, purity .95-.97, coverage .99, readouts = true-mask fits; qtoyflex beats oracle masks; nulls 0 sites); traps: noise-fitting, support-dof prior, frozen-lambda death spiral, union merging"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-07T22:02:47.635Z
---

Side session sparse-ppca, 2026-11-07. Handoff: `hetero_tree_20261028/handoff/sparse-ppca/{VERDICT,STATUS}.md` + scripts/; runs in
`math_reworks/sparse_ppca/runs/` (x10fix_p, x10null_p2, x5fix_p2, x5null_p2, qflex_p2, qnull_p3; comparators x10fix_ta*b).

Python prototype `sppca.py` on the latent framework (lp 20 / box 70 IgM, lp 15 / box 36 qtoy), particles x Fourier-cell sparse
operators + numba kernels, 12 nuisance templates (rigid tangents etc.) as latents.

- **IgM 10 % fix:** 5/5 arms, one 2-mode block each. Purity .95-.97, coverage .99. R2 / LDA / GMM .39-.52 / .29-.39 / .28-.39, equal
  to the true-mask focused fits; GMM above T1 on every arm.
- **IgM 10 % null:** the continuously moving arms are found, with no discreteness (held-out 2- vs 1-Gaussian gain < 0).
- **IgM 5 %:** 4/5 arms at T1-level readouts, 1 weak junk block. The 5 % null gives 0 sites.
- **qtoyflex:** closed AUROC .897 / holo .926, against oracle-mask fits .852 / .851 and flex's node latent .757 / .861. qtoynull
  gives 0 sites.

**Why:** BBP. An arm's whitened eigenvalue (~2.3) is below a global PCA's noise edge (~4.8) but above an arm-sized support's (~1.0),
so global latents and rotations cannot see arms; local seeds plus sparsity can.

**How to apply / traps (each measured):**
- Unregularised ML PPCA fits noise.
- The per-shell prior must count only the mode's support dof, else localised modes ARD-collapse.
- A group lasso with a frozen lambda spirals W to 0. Use select-then-refit with lambda recalibrated each iteration.
- All seeds at once merge into union supports. Use greedy pursuit with held-out gains vs solvent controls.
- The rotation step needs one lambda per group, else blocks merge.
- Keep a leave-one-block-out check, else the last block of a null survives.

**Knob-free stop (17:00, STATUS.md "Coordinator ask 1 / 2", runs `*_p3`):** pursue while the best seed beats 2 SE and all 12
solvent controls; keep blocks by the same test left one out; sites = blocks joined by support overlap above hypergeometric chance
+ 2 SE; states by a TIED-covariance GMM (a full GMM on a multi-block site splits heavy tails: a2 .003 vs tied .290). Stops by itself:
10 % fix 8 blocks -> 5 arm sites, 5 % fix 4 arms + 1 weak, qtoyflex 1, 5 % null / qtoynull 0. The nulls' seeds DO beat the solvent
controls (5-37 nats vs <= 0.9); the support-dof shell prior drains their grown blocks to 0, so LOBO drops them. The 2 % cut's
"sixth candidate" was a second seed on arm a5. The 2- vs 1-Gaussian gain is a readout only (Part 4: as a calibrated test it
detects non-Gaussianity). 10 % null knob-free rerun pending (13 blocks by round 13).
**Next (user, 2026-11-07):** "Port sparse PPCA first" (native opt-in flex_pca path, knob-free stop replacing the 2 % cut; given to
the sparse-ppca2 session), then the end-to-end chain. Its sparse global coordinate goes to real data "After the port": r131 roots
(both seeds) and the 692k root with |CTF| ([[rnasep-denoised-stacks-are-phase-flipped]]), against the seam-free .69-.72.
Lever 3 on these supports: 10 % 2-state arms sit at the masked rounds' fixed point (~.35-.37); 5 % fix product .104 at round 1,
above the true-mask chain's .098 at round 3.

Related: [[site-dictionary-global-basis-verdict-2026-11-07]], [[igm-focused-flex-mask-is-bottleneck-2026-11-07]],
[[site-mask-rework-direction-2026-11-07]], [[flex-pca-sparse-l1-mstep-verdict]].
