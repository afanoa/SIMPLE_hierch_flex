---
name: flex-delivered-state-maps-one-process-2026-10-24
description: "flex_pca's delivered state maps (deliver_state_maps -> reconstruct3D without nparts) run in ONE process: 34-47 % of a depth>0 fflex step (1.06 ms/particle at 8 thr); distributed replay gives the same maps (corr 1.000000000), r1 599 -> 189 s; the root's rec covers only updatecnt>0 rows (rnasep.simple: exactly 25,000), not a cap"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-24T05:39:33.180Z
---

2026-10-24, flex / fflex speed-up profile for the coordinator (read-only + replays in batch2 `rnasep_klem/fast_tree/prof/`:
drec_r1, drec_cahra, fshm_r22; scripts fast_drec_v1.sbatch, fast_fshm_v1.sbatch).

- `deliver_state_maps` (simple_flex_pca_project_gateway.f90) builds a reconstruct3D cline with no nparts, so the K-state maps are
  one in-process job even in a distributed flex run. Below the root (oneflex: fflex only) it reconstructs every node row:
  ~10 s + 1.06 ms/particle at 8 threads (m33: r1 482k 521 s of a 1,115 s fflex; r11 366 s; r2 258 s).
- Replay with nparts on the step's own binary and published project: r1 (12 parts) 599 -> 189 s incl. volassemble; CAHRA root
  fflex (30k C11, 4 parts) 467 -> 195 s. Maps equal to the step's: corr 1.000000000, max rel diff <= 2.6e-5 (summation order).
- A second distributed r1 run was byte-identical to the first (32 of 32 files). Patch for the user (not applied):
  `.claude/patches/flex_deliver_maps_distributed.patch` on 056895a5c.
- Validated in place with TEST build batch2 `~/SIMPLE_fastl1` (056895a5c + patch, a plain source tree from git archive, built
  by fast_tree/fast_l1build_v1.sbatch in 4 min on 16 cores; fast gate 15/15): through flex_pca, maps = old in-process rebuild of
  the same weights (<= 2.5e-5); r1 maps 703 -> 156 s, CAHRA full fflex step 720 -> 470 s.
- In seeded tree runs (SIMPLE_SEED exported) the root flex and fflex draw the same sigma2 bootstrap sample: m33's two
  sigma2_state.bin are byte-identical, so fflex reusing flex's sigma2 is identical.
- Pre-launch check (job script fast_tree/fast_updcheck_v1.sbatch): `simple_private_exec prg=print_project_vals projfile=<p>
  oritype=ptcl3D keys=updatecnt | awk '$2>0{a++; if($3>0)u++} END{print "active="a+0, "updated="u+0}'`.
- The root's rec (and the root mean) reconstruct only rows with updatecnt > 0 (sample4rec): rnasep.simple has exactly 25,000,
  so plan v2's "flat 25k root rec" is that subset. A 4.5M input with all or no rows updated makes the root flex and fflex each
  rebuild 4.5M rows in one process (~80 min each) unless the rec is distributed.
- Part rounds re-read and re-prep every particle (RESIDENT PLANES OFF without cache=yes): one pass at r22 (40k, 8 thr) is read 21 s
  + prep 4.7 + project/solve 1.8 + insert 15.4. cache=yes as built is master-built (the master reads every particle serially), so
  it is a loss in the distributed tree.
- The flex master runs at the step's nthr (8 in m33), the same as a part; DECONV and the rec scale ~1.9x for 2x threads (CAHRA p4
  vs p8).
- Determinism leaks found: schedule(dynamic) into per-thread sums (E-step loop, em_calibrate_noise_prior), reduction(+:) on reals
  (XD fit loglik that also drives its stop, deconv noise scale, crossfsc), the sigma2 bootstrap sample from /dev/urandom. The
  six code sites are exactly the port's flex_determinism_on_056895a5c.patch; the bootstrap sample needs SIMPLE_SEED.

**Why:** the user's 12 h target for ~4.5M; flex_pca is the user's code, so diffs go to the user first.
**How to apply:** for a flex wall-time question, look at the delivered-map rec first; check the input's updatecnt before a big
root. Related: [[hetero-tree-4p5m-12h-serial-profile-2026-10-23]].
