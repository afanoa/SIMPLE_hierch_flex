---
name: klem-state-graph-umap-and-real-halves-2026-10-19
description: "state_graph.py (cross-half similarity + shared-particle edges) and umap_states.py (UMAP of cross-half memberships: a starfish, arms = states, mixed core; the old island UMAPs were built from per-board responsibilities); real RNase P gold-standard halves rnasep692 (eo split, |CTF| geometry check FSC 0.68 at 10 A), dtree gs1 running"
metadata:
  type: project
---

2026-10-19 (user: build the graph view; "can we still expect a umap to look like this" = the old island UMAP).
- toy55/state_graph.py + state_graph.sbatch (MODE TAG FIN): nodes = final states, blue edges = cross-half similarity of
  state deviations (39-10 A), red dashed = shared particles (cosine overlap of soft memberships from one cross-half
  E-step), MDS layout, after-the-fact average-linkage dendrogram; saves memberships WA/WB in the npz.
  Six-conformer 700k: the two 42 states link to each other and to 46; the mixed state is a hub sharing particles
  with 49, 46, 47 and 29+27. 55_200k: junk states a tight cluster, real states a chain (continuum), one tiny state
  not reproducible.
- umap_states.py + umap_states.sbatch: UMAP of sqrt(memberships). Six-conformer: one connected starfish, an arm per
  state from a mixed core, truth enriched per arm but mixed; NO separated islands. The earlier island UMAPs
  (global_embedding.py, umap_leaves.py) concatenated each particle's responsibilities over the boards it passed
  through, so different branches sat on disjoint coordinates: islands by construction.
- Real halves (real_halves.sbatch): flags==1 rows split by meta eo (345,999 / 346,393), stacks written row-ordered
  through stkind/indstk, projects imported ctf=flip (kernel snapshot build), run dirs synth_rnasep692_A/B, caches
  benchcache_rnasep692. geometry_check.py (40k per half, bench insertion): |CTF| FSC 0.97/0.91/0.86/0.68 at
  19.5/14.6/11.7/10.2 A; signed CTF 0.80 at 14.6 then 0.10 -> real data needs --ctfabs 1 / CTFABS=1.
  TRAP: indexing an NpzFile per row re-reads the array each time (first prep crawled): load once into a dict.
- Real tree: dtree_v12 MODE=rnasep692 TAG=gs1 BOARDJOB=dtree_board4 CTFABS=1 PSCALE=1 BENCHJOB=splitbench_node5
  BENCHCTF=1 NOGT=1 (job 54008442); then tree_finalize3 (CTFABS=1 NOGT=1), state graph EXTRA="--ctfabs 1" NOGT=1,
  state_recs.sbatch (CLOSED=maps160/docked_all/conf18_emd74209_m160.mrc, CONSMAP=rnasep_states/25_refine3D_states/
  rec_final_state01.mrc) for the closed-arm readout. Related: [[klem-six-conformer-700k-and-junk-flag-2026-10-19]].

REAL ROOT (gs1, 22:4x): K* 4 = r.1 163k/148k, r.2 77k/92k, r.3 66k/65k, r.4 junk 40k/41k (free share 0.97/0.98, flagged).
Bench peak MaxRSS 114 GB of 120 at 346k per half (loader keeps every temporary) -> final pass / graph jobs on the real
set at --mem=200G (nodes 256 GB). Closed proxy = f692c fam1 zoomregionmotion_1_s05 sub-states 1+3 (0.92 closed-arm,
5.2k/5.6k per half, 1.5 %): ~60 % of them go to r.2 in both halves (3.9 % there), junk leaf holds none.
Arm-focused zoom on r.2 per half (arm_zoom.sbatch H= SEL=node:r.2; mask dtree_rnasep692_gs1/arm_mask.mrc from
make_arm_mask.py, 3.3 % of the box), jobs 54008984/85; report via zoom_report.sbatch SELN=node_r_2.
ARM-MASK ZOOM ON r.2 = NEGATIVE (22:59): one supported sub-state per half (closed-arm 0.51 A / 0.73 B, known-closed share
1.9 % / 7.9 %); closed_arm_cells.py gives the SAME occupancy for every cell (0.59 A, 0.62 B): the arm-travel-mask board
finds no arm structure (as the 2026-10-17 arm boards). Next = the route that worked before: region_route.sbatch (H, TAG,
NODE): derive_regions2 on the tree's own r.2 board (klem_fortran2npz first: dtree boards keep only .bin), region with the
largest arm-mask overlap as focus, region board -> node_states2 level 1 -> rec + arm_readout + proxy share -> zoom on the
most closed level-1 state (7x7 + 1) -> sub-states -> readouts. Jobs 54009008/09.
REGION ROUTE ON r.2 = NEGATIVE TOO (23:20): the tree's r.2 board kept only 9 (A) / 7 (B) reproducible cells; derive_regions2
regions OK (region_1 covers 0.74 / 0.43 of the arm mask) but the region board -> ONE supported level-1 state per half
(closed-arm 0.51 / 0.63) and the zoom -> ONE sub-state (0.53 / 0.62), known-closed share 2-5 % (no enrichment).
Reading: f692e's boards worked at FAMILY scale (~300k, both halves, refined poses); a gold-standard half-branch has
~80-90k -> node_states2 --support finds one state. Recommendation given to the user (2026-10-19 23:3x): discover on all
particles (f692e route), VALIDATE with the cross-half final pass on the independent halves (f692e leaves + closed
sub-states as candidates via meta_rows), instead of discovering inside half-size trees.
