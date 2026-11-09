---
name: flex-basis-and-z-signs-arbitrary-across-runs
description: "TRAP: flex's written basis volumes (merged/pc/polished/fitB_pc*.mrc) and the embedding's z carry independent arbitrary per-component signs across runs/seeds; align by per-particle agreement of W z, not by state maps"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-09T05:29:20.205Z
---

Found in xd-repro (2026-11-09) on the r1 node-levers replays. Two seeds with the same latent (identical XD fits) differ in sign,
independently, in two places:
- the merged basis flips components 4, 6, 9 and 10 (same-index corr exactly -1);
- z flips components 7, 8 and 9.

The pc, polished and fitB files flip different components again; pcit000_ flips none.

Two consequences:
- Unsigned deviation maps W(mu_k - mu_bar) correlate at only .32-.39 between identical fits.
- Reading the signs from flex's state maps fails. Projecting recvol_stateNN.mrc on the basis against flex_pca_state_targets.txt
  gives |r| .03-.98 and is wrong at 3 of 10 components.

**Why:** any cross-run comparison of latent-model maps, bases or latent coordinates silently mixes sign conventions.

**How to apply:**
- Set each run's signs against a reference run, by per-particle agreement over the shared rows:
  s_j = sign(sum_k (W_B' W_R)_jk (Z_B' Z_R)_jk).
  At r1, |cos| was .92-1.00 and the merged basis is orthonormal.
- Code: `math_reworks/xd_repro/scripts/xr_repro3.py`.
- The median |corr| of z across runs hides the flips, so don't use it as a frame check.
