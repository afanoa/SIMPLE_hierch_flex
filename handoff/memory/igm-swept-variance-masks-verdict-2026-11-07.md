---
name: igm-swept-variance-masks-verdict-2026-11-07
description: "VERDICT guide 7.8 levers 1+7 (IgM fixtures): swept-variance masks from latent-binned reconstructions reach true-mask numbers at 10 % where islands find the arm, fail at 5 % (bins carry the flex modes' selection bias); lower per-site lp never helps, keep lp 15"
metadata:
  node_type: memory
  type: project
  originSessionId: 90c3ae86-904a-4714-a0d5-877dd30cbfc6
  modified: 2026-11-07T19:36:53.650Z
---

2026-11-07, side session of the site/mask coordinator. Handoff: `hetero_tree_20261028/handoff/igm-masks/VERDICT.md`
(tables `results.md`). Runs and masks: `igm_cryolo/tree/fixture_igm{,10}/igm_masks/`.

**Lever 1.** Path: trimmed island fit -> 6 z1 quantile bins -> reconstruct3D halves at the root's poses -> cross-half bin
variance -> mask -> refit.
- 10 %, island fix arms, rule vs > 100 x solvent sigma (anchor-connected, data-only): a1 R2 .547 / LDA .433 vs true .556 /
  .444; a3 .474 / .341 vs .455 / .320.
- 10 % null: angle R2 .277 vs true .274; tied gain falls to the true-mask level.
- Iteration 2 rescues a weak transferred start when the signal is strong (10 % a2 .197 -> .419), not a4.
- 5 %: no variant beats the trim (only a1 gains).

**TRAP.** Binning by a particle's own focused coordinate bins its noise along the flex modes. Every bin map then carries the
mode shape inside the start mask, both halves alike, so the cross-half variance keeps it (5 %: 1400-2400 sigma inside vs
<= 196 outside). A peak-relative threshold (0.2 x max) re-makes the start mask, like T1b's `self`. A band split cannot
debias here (no state signal at 11-7 A). Also: z1 bins build the wrong mask when the start fit's z1 is not the arm (5 % a2,
a4).

**Lever 7.** Mean fix R2 falls monotonically lp 15 -> 30 at both scales and the null is flat. Lower bands heavy-tail the
null coordinates and erase the 10 % fix-vs-null tied-gain separation. Keep lp 15.

Other findings:
- The 10 % fix root's recB0 islands cover a1/a3 only: discovery is still the gate.
- The 3-D GMM NMI is fragile (.000 vs T1's .300 on coordinates equal to 3 decimals).
- A recipe with the flex parts run in-job (qsys local, SIMPLE_SEED=1) reproduces T1.

**Why:** answers whether a data-only mask can replace the true arm mask for focused flex (T1's bottleneck,
[[igm-focused-flex-mask-is-bottleneck-2026-11-07]]).

**How to apply:**
- Port candidate: the swept-variance refinement after the trimmed island, at real-Fab SNR only, with a noise-calibrated
  threshold.
- No truth-free gate for when to adopt it was found (outside-support variance overlaps between gains and losses).
- Don't propose lower per-site bands again.

Related: [[site-mask-rework-direction-2026-11-07]], [[igm10-fixture-real-fab-proportions-2026-11-07]].
