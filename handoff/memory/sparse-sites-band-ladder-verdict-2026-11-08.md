---
name: sparse-sites-band-ladder-verdict-2026-11-08
description: "BUILT: sparse_sites band ladder (pursued-gain argmax over quarter-octave caps, one octave down), patch band3 md5 b5d9c457 on wip 7ae8ba4b1, build wip_builds/band3_7ae8ba4 (gate 15/15); IgM 10 % fix 18.2 A 5/5 arms R^2 within .002, null 18.2 A (flat plateau, no disc), 5 % fix 25.2 A, qtoy/r131 full band (noise model on too); TRAPS: a cap must be a fresh run (raw samples + RNG reset); dropped rounds drift the kept blocks; cost 1.6-2.3x per site node, 4-4.7x null/rising"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-08T15:20:55.494Z
---

Side session sparse-ppca, 2026-11-08. The coordinator's band rule built into flex_pca sparse_sites
([[sparse-sites-native-port-2026-11-07]], [[band-must-be-data-derived]]). Handoff `hetero_tree_20261028/handoff/sparse-ppca/band_ladder.md`
(also the shared notes folder as `sparse_sites_band_ladder.md`). Runs `math_reworks/sparse_ppca/native/lv3_*`.

**Status:**
- Delivered: `sparse_sites_band3.patch`, md5 b5d9c457, on fork wip `7ae8ba4b1`, which carries the noise model `03506060d`.
- Build: `wip_builds/band3_7ae8ba4`, fast gate 15/15. It needs `--constraint=x6342`.
- Nothing committed; the coordinator commits after its gate.
- v1 and v2 builds (`band_5b17de8`, `band_7ae8ba4`, `band2_7ae8ba4`) are superseded.

**Results (default whitening):**
- IgM 10 % fix: 18.2 A (ring 36; 17.7 A is not a ladder point), 8 blocks in 5 sites, one per arm, every arm's R^2 within .002 of
  the scan's best.
- IgM 10 % null: 18.2 A, not the scan's 20.4 A. The pursued gain is flat over 17.7-21.8 A. No site is discrete.
- IgM 5 % fix: 25.2 A (ring 26), 5 arms. Same R^2 as the lp-20 lattice's 20.4 A runs within .003, a2 +.016.
- qtoyflex and r131 pick the full band. qtoynull gives 0 sites at all 5 caps.
- r131 with `SIMPLE_FLEX_NOISE_VAR=yes` also picks the full band (9 blocks, closed best block .681, span .765).

**Traps:**
- **A cap must equal a fresh run at that band.** Restricting the full operator's whitened samples, or letting the RNG stream run on
  across caps, gives a different and weaker pursuit. On IgM ring 36 ran 19+ rounds of 1-2 nats.
  - v3 keeps the plane pass's raw samples. Each cap re-assembles them with its own beta / sig2 / amplitude (`ops_band`) and resets
    the RNG (`get_rnd_state` / `put_rnd_state`).
  - qflex's ring 13 then equals the scan's `kqflex_c175` to every digit.
- **A dropped round is not undone.** GROW_ITS refits every block with the candidate in, and `drop_block` removes only the
  candidate. Many drops in a row drift the kept blocks: r131 noise-on ring 13 fell from .212 to .015 LL/ptcl over 11 drops. This
  predates the ladder and is flagged as a follow-up.

**Cost:**
- The ladder costs 1.6-2.3x the node-band run per site node and 4.0-4.7x per null or still-rising node.
- On IgM 45-57 % of the time goes to caps that keep nothing. Their seeds pass round 1 and the grown blocks decay to zero, so a
  round-1 screen saves nothing there. It saves 73 % on box-36 null nodes.
- Caps are independent, and caps 1-2 always run, so they could run as two jobs at once.
- 692k root extrapolates from 130 min to 3.5-4.3 h. Peak memory is 1.3-1.55x the single-band peak.

**How to apply:**
- The site chain's IgM trees run on band3 plus the v8 candidate.
- Judge the null fixture on discreteness, not on its exact band.
- Re-read r131's pick if the noise model becomes the default.
- 2026-11-08: commit HELD. The user asked for the band inside the model instead of the ladder ([[band-must-be-data-derived]]).
  - sparse-ppca2 prototypes it against the ladder's seven runs. The ladder is committed only if the in-model band fails.
  - The site chain's b11 (b10 + band3) IgM trees run meanwhile.
- 2026-11-08 later: the in-model band failed on IgM in four forms ([[sparse-sites-inband-band-verdict-2026-11-08]]); the ladder stands.
- 2026-11-08 evening: the LADDER IS FINAL. Commit f4be24169, branch sparse-sites-band-ladder off wip 20e972dfd (worktree
  /Users/atifao/SIMPLE_ladder): band3's rule, plus two cost changes.
  - Cap 2 runs as a queued sub-job beside cap 1.
  - A round-1 screen for caps 3+ applies only when the node band finds no seed.
  - Gates: build + fast gate (55086887), then the 7 ladder fixtures must reproduce every pick, block count and gain. Only
    then cherry-pick to wip.
  - TRAP: cap sub-jobs need x6342 (USE_ARCHOPT). The parent job sets SBATCH_CONSTRAINT=x6342, and tree launch scripts must
    set it too.
- TRAP (IgM check, 2026-11-08 evening): eb1d99ede's cap-2 takeover renames a folder, and BeeGFS refuses the rename. The winning cap's
  sites were lost (0 sites on IgM 10 %). Fixed in e4fbef54e (on wip cc4af8a10), now re-validating. Never push eb1d99ede.
