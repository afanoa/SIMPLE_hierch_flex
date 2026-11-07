---
name: flex-pca-fork-map-2026-09-15
description: "the three flex_pca lines since 0fcd70404 and where each is built; flex-delivery (old layout, fixed states, no PCG) vs flex-pca-pcg (new layout, PCG, hard k-means states); DIRECTION DECIDED 09-16: port PCG onto flex-delivery, file consolidation deferred"
metadata: 
  node_type: memory
  type: project
  originSessionId: ab9634b6-edc2-4f26-be6e-dc10b4a08ba1
  modified: 2026-09-15T15:50:43.510Z
---

Fork point 0fcd70404 (2026-09-09). Three lines, none rejoined as of 2026-09-15:

- `flex-delivery` a94d6c790 (33 commits): OLD 33-file layout (weights/gmm/targets, em_polar...).
  The 09-12 cleanup tranche (d234f83df..a94d6c790): validated soft-weight placement, occupancy
  floor, 84->14 env reads, `rec_states` CLI key. Verified bit-for-bit vs 10028 `bare_b64_lp16`.
  NO PCG (rec_backend ignored by flex_pca there). Built: box `~/wt_res/build`, batch2
  `~/SIMPLE_delivery/build` (09-15, gcc 15.2, ~5 min at -j24).
- `flex-pca-pcg` 8e081f9d9 (96 commits): NEW 26-file layout = the Mac consolidation snapshot
  8f8f5ba5a + rec_backend=pcg + pcg_mskfile + smpd_target + 27x speed work + master merged to
  d8ae11d15. Still 79 env reads and HARD unit-weight k-means states (neff==particles). This is what
  every batch2 sp_*/cn_* script runs (`~/SIMPLE_pcg_opt10..12`, `~/SIMPLE_pcg_mem`).
- `flex-pca-merge` 275844acb: master port WIP (Slice A only); `pcg_mskfile` still bound 'B' so any
  MRC mask is refused at the cmdline check. Retire.

The Mac's dirty `/Users/atifao/SIMPLE` tree == 8f8f5ba5a for src/main/flex+commander+params+ui+tests
(diff empty). batch2 `~/SIMPLE_pcg_mem` = 8e081f9d9 + UNCOMMITTED "band boxing" part-format change
(+437/-72 in em_estep/rounds), built 09-15 00:05.

**Why:** the cleanup was done on the old layout because the consolidation was never committed on the
box; merge-tree gives 10 content + 4 modify/delete conflicts, so `git merge` is not the tool.
**How to apply:** reconcile by re-applying the delivery tranche knob-by-knob ON flex-pca-pcg (branch
`flex-pca-unified`), state stage first (d234f83df, d0acf7f58, 1889ff82e). Plan + verification matrix:
`doc/refactoring_notes/flex_pca_branch_reconciliation_2026_09_15.md`; harness `scripts/flex_pca_verify/`
(= batch2 `~/flex_verify/`, box `~/verify_delivery_20260915.sh`). Run roots
`.../10028_simple/flex_verify_20260915/<arm>` and `.../10180/flex_verify_20260915/<arm>`.
See [[flex-pca-default-recipe-verified]], [[flex-pca-pcg-branch-state]], [[flex-pca-complement-block-and-step0]].

**User direction 2026-09-15:** the spliceosome recipe (paired probe-only, pcg+automask basis, box 64 lp 20,
46k selection) "works mostly as expected"; PCG is a non-issue. The open problem is STATE PLACEMENT and it
must be AUTONOMOUS: decide state_axis (0 vs an axis), how many PCs (nkern), without expensive
reconstructions (decide in the latent on the embedding cache; reconstruct once). Backend split: pcg for the
basis, gridding for the final state maps (keep all 4 combos as options later, not implemented yet).
Verification relaunched 13:24 (batch2 jobs 53213368-79 + box chain) after the BeeGFS outage.

**2026-09-16:** the box's `~/SIMPLE_fast` (NOT a git checkout, built 09-15 23:51, runs trpm4) = 8e081f9d9 +
band-boxed part files (formats 10/11) + typed `rec_states=yes|no` and `rec_states_backend=gridding|pcg`
(default gridding for the state maps even under rec_backend=pcg). It is the user's current best and is
unversioned: commit it onto flex-pca-pcg FIRST (diff saved in the session scratchpad as
simple_fast_vs_wtpcg.diff, 1126 lines). Implementation plan for the unified branch:
`doc/implementation_notes/flex_pca_unified_implementation_plan_2026_09_16.md` (defaults: pcg+mask basis,
gridding states, smpd_target primary with a cost guard, warm-scaled PCG start, autonomous placement in the latent).

**Direction reversed 2026-09-16 (user):** the 32->22 file consolidation (8f8f5ba5a) is pure motion and is to be
saved for a later, separate mechanical commit. So the unified branch is based on flex-delivery (old layout,
knobs already gone) and the PCG line + SIMPLE_fast delta are ported ONTO it (additive: pcg.f90, rec3D_pcg.f90,
hooks in shared files). Order in the plan doc section 2: master sync to d8ae11d15 first, then PCG, part
formats, solver speed, smpd_target, rec_states_backend, placement, speed, master tip, consolidation last.

**Defect found 2026-09-16:** flex-delivery and flex-pca-merge write the boosted master thread count (nparts*nthr)
into the part scripts (`--cpus-per-task=160` at nthr=16) -> sbatch refuses, master dies after 5 retries. Cap the
master with SIMPLE_COV_MASTER_NTHR where it exists (pcg) or keep nparts*nthr <= 96 on csbdevel; the unified branch
must derive the part job_descr nthr from the worker nthr.

**User decisions 2026-09-16 (placement):** (1) the per-state eo-FSC low-pass on delivered state maps is ALWAYS on
(no STATE_EOFILT switch; it supersedes the 09-07 "consensus filter looks right" note). (2) On the 10180 46k
spliceosome the default state_axis=0 (deconvolved-mixture macro-clusters + k-means regions in the standardised
latent, S_P_grid: 3936 / 4 x ~10.7k, three states filtered to 25-29 A) does NOT capture the conformational
states; equal-mass partitions along PC1 do. The autonomous rule must pick PC1 equal-mass there from the latent
alone, with no reconstructions. See plan doc section 3.
**User decision 2026-09-16 (outputs):** deliver only the ranked principal components + the filtered state volumes
(with eo halves); initial bases (pcit000), per-fit fitA/fitB and merged sets become diagnostics-only (plan doc 7a).
**User decision 2026-09-16:** n_probe_iters default = 4 (the commander currently sets 5 when the key is absent; params field is 0).
**User decision 2026-09-16:** sigma_est=global must not be a required knob on selected subsets; flex_pca auto-falls back
to the pooled spectrum when a sigma2 group is empty (plan doc section 7); shared-reduce fix is a Hans conversation.
**User decision 2026-09-16:** box_crop + lp fold into ONE key: `lp` (resolution of the variance); smpd_target=lp/2.5 and box_crop derived (box_crop expert override only).

**Reconciliation started 2026-09-16 evening.** Frozen results tree: tag `flex-pca-verified-2026-09-16` = 8e081f9d9
(batch2 `~/SIMPLE_pcg_opt10` build; NEVER rebuild/checkout there). flex-pca-pcg now = 3e80b4265 = tag
`flex-pca-pcg-final` (band boxing 0dae9c8dd, rec_states keys dac07b5c1, figures+one-file cache 3e80b4265);
gate build in batch2 `~/SIMPLE_pcgfinal`. Unified branch `flex-pca-unified` = flex-delivery a94d6c790, tag
`flex-pca-unified-base`, Mac worktree `/Users/atifao/SIMPLE_unified`, pushed to batch2 and box. Next: port
step 1 (master sync to d8ae11d15 + sigma fallback) per the v2 plan section 3.
