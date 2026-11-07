---
name: hetero-tree-schedule-cuts-fail-unconverged-drift-2026-10-22
description: "TRAP (fast-tree, r1): hetero_tree's stage A+B partition is NOT converged (10.8k -> 38.1k over 6 its); cheaper schedules fail a v47 LABEL standard by construction. RESOLVED 2026-10-23 by scoring against TRUTH: 1+1 per stage wins on both toys, 3+3 on CAHRA by 0.03; the user made 1+1 the default"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-23T18:57:28.261Z
---

2026-10-22, fast-tree session (runs `rnasep_klem/fast_tree/shift`, `/site`), all r1 replays from v47a's inputs. The noise standard
is v47a vs v47b: site Dice 0.997, seeds NMI 0.980, stage-B NMI 0.948-0.958.

Minority size at the r1 root:
- seeds 10.8k
- stage A its 1-3: 24.6k, 31.1k, 34.3k
- stage B its 1-3: 35.2k, 37.1k, 38.1k (state overlap 0.98 at the last)
- refine=shc after 12 its: 40.0k (overlap 0.995)

Cheaper schedules tried (each fails the standard):
- stA0 at 1 iteration: discovery picks other sites (Dice 0.87).
- stA0 at 2 iterations: same sites, Dice 0.987, but seeds NMI 0.94.
- Stage B straight from the seeds: NMI 0.62.
- refine=shc: NMI 0.55, and 2.8x slower than prob_neigh.
- trs=0: NMI 0.84.
- npeaks_inpl=2 in every step: NMI 0.91.

What passes:
- npeaks_inpl=2 at stA0 only. stA0 feeds only site discovery; site_rank, the seeds and stage A all start from the node's input poses.
- The calc_athres threshold cache, identical by construction: `.claude/patches/threshold_cache_on_51c8cdac4.patch`, -19 % at r1.
- site_rank is deterministic: my build on v47a's mask gives seeds NMI 1.000.

**Convergence test, done 2026-10-23** (`fast_tree/conv`, merged build, vm8a's seeds; every iteration forced with minits):
- The partition converges slowly: overlap >= 0.995 by stage-B it 7-10, minority ~41k at 15 its, still +50-110 per iteration.
- Same schedule, two runs: labels NMI 0.93. Only 0.6-0.85 % of labels differ, all in the lowest ~5 % of state margins.
- Different schedules converge to different labels (NMI 0.49-0.72), but the maps are within noise of each other.
- So the 15.6 % "boundary" is schedule dependence, not replicate noise.

**Margins** (from the last prob tables, dist_part_neigh_*.dat):
- Low margins predict flips: AUC 0.98 within a schedule, 0.84-0.91 across schedules.
- The last-iteration changes (mi_state) do NOT calibrate replicate flips: 5x too few at convergence, 4x too many at 3+3.
- The objective's own scale: LLR = margin x wsqsums / polar oversampling.
- The user approved margin-weighted children (2026-10-23). The recipe went to the port.

**Why:** the user wanted 2-3 h for the 692k tree from a cheaper refine schedule. The reference endpoint is wherever the 6
reassignment iterations stop, so the standard is not schedule-invariant.

**Resolved 2026-10-23 by scoring against truth** (`fast_tree/tset*`, design note "Contrast-on rerun"). With the per-particle
contrast on (m8c27), NMI against truth, (a) 3+3 vs (b1) 1+1 per stage:
- 2-pole toy: 0.275 vs 0.461.
- 3-pole toy: 0.284 vs 0.485.
- CAHRA c11: 0.363 vs 0.329. The replicate pair reproduced both to < 0.001 NMI (labels agree > 99.9 %), so this gap is real.
- b1 costs a third.
- The user chose ONE PASS PER STAGE as the hetero_tree default; the port builds it.

**How to apply:** don't judge a cut schedule against v47a/v47b labels. Score schedules against truth (toys, CAHRA) with
replicates, as above. Re-opening the schedule needs new truth evidence, not label drift. Per-node time is dominated by fixed per-iteration cost
(~60 %) and the root by 8-thread site_rank / state_split. Related: [[fixed-pose-node-mask-em-is-not-refine3d-2026-10-22]].
