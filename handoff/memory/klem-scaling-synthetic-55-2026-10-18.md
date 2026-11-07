---
name: klem-scaling-synthetic-55-2026-10-18
description: scaling test of the state count vs particle number: the full 55-conformer RNase P landscape simulated at 200k / 700k / 4M (run_toy55_scale.sbatch, MODE=55_<tag>, poses with replacement via make_oris55r.py) through the FULL chain (synth_tree_prep_full.sbatch), readouts = chain landscape counts, synth_truth purity, readout55 strict recovery
metadata:
  type: project
---

User 2026-10-18: "how do we avoid splitting so much that resolution suffers ... if I scale to 4M RNase P particles am I
going to get the same number of states?" Answer: states exist only by support (cross-half reproducible, separable,
12 A floor, >= 1000 ptcls) + the data-derived merge; the count should rise with N and stop at the data's conformers, not
at a tree-set number - unmeasured until this test. Setup: all 55 docked apo conformers (maps160/docked_all) at the
paper's image-count proportions, real poses + CTF (meta_all refined rows; with replacement when N > 692k), calibrated
bfac/snr (ah calibration), 10 % blob + 10 % broken junk, refine3D prob lp 8 from 10-deg perturbed poses; sims
sim55_200k / sim55_700k / sim55_4M (jobs 53988511-3; 4M = ~410 GB images, beegfs has 4.2 PB free), projects
proj55_<tag>, pose tables lat_55_<tag>. Chain = f692_chain9 production presets (no DEV subsampling: DEV would cap
at 150k and void the test). Readouts: chain "similarity 0.5-0.9: N states / identified", synth_truth.sbatch (per-leaf
majority conformer + purity: the right metric for discrete conformers), readout55.sbatch (strict recovery of the 55).
Expected: states vs N curve against the 55 line; particles per state; conformers per state. Watcher chains the stages.

**Over-splitting tests (user: "how are we sure this method wont oversplit small data sets", 2026-10-18 ~20:45):**
null_50k / null_200k = ONE conformer (49) + 10 % blob + 10 % broken junk (null_conformers.csv), two_50k = conformers
18 (closed) + 49 at the paper's weights (two_conformers.csv); run_toy55_scale.sbatch now takes CSV and simulates
only the csv's conformers (jobs 53992535-7), then synth_tree_prep_full -> full chain -> synth_truth. Right answers: 1
state (+ junk) at both sizes; exactly 2 states with high purity. Every extra supported state = over-split; the truth
readout says what it is made of (view artefact / junk / split of one conformer).

**200k point (synth_55_200k, full chain 20:11-21:33 = 1.4 h):** 2 families (25 + 24 states), cell-states 15 / 24 / 71 /
215 at similarity 0.6-0.9 (real 692k: 75 / 109 / 156 / 253), docked 10: 7-9 identified; truth: leaf purity 0.10 (leaves
are mixtures of many conformers at ~3.6k particles per conformer), 7 of 100 labels majority of some leaf; strict
recovery 8 / 55 (mean best-match corr 0.64). 700k chain started 22:49; 4M simulated 22:13 (merging/refining); null_50k
and two_50k chains running since 22:13; null_200k in its refinement.

**Diagnosis (user, 22:5x: "It seems like its oversplitting. Its hard to recover enough signal if states are too
small"):** at 200k the leaves pass support (reproducible, separable) with purity 0.10 = reproducible MIXTURE slices
(separability is not purity). Fix under test: merge_data_derived.py `--floor auto` (env FLOOR=auto, MARGIN): the
difference of two states must be resolved at the POOLED state's own FSC-0.143 resolution x margin (band 1.25 d .. d),
not at a fixed 12 A; mixture slices differ only at low resolution -> merge; closed/open differ at the map resolution ->
stay. consolidated_purity.py reads merged states against the truth; merge_rule_test.sbatch compares floor 12 vs auto on
synth_55_200k (job 53994044). pooled_dd.sbatch maps FLOOR=auto to CUT=12 for the directory name
(consolidated_12_rule_auto).
Merge-rule test on synth_55_200k (13 leaves): floor 12 -> 10 states (median 5.9k ptcls), purity 0.07, strict 3/55,
docked 6/10; floor auto -> 4 states (median 24k), purity 0.06, strict 2/55, docked 3/10. Purity cannot decide here:
55 conformers at ~3.6k each cannot be pure in 4-10 states (ideal majority share ~0.2). Right readout for a dense
landscape = coherence: do the conformers co-occupying a state form a contiguous region of the landscape (mean pairwise
docked-map deviation similarity vs random sets, state_coherence.py / coherence.sbatch)?

**Over-splitting tests, leaf level (2026-10-18 23:xx):** null_50k (one conformer + junk): every region board yields ONE
supported leaf (one board two), purity 1.00 - the support test does not split homogeneous data. two_50k (18 + 49): leaf
purity 0.96, both conformers the majority of leaves. The chain's "similarity 0.7: 127 states" line counts raw lattice
cells, not supported leaves - never quote it as a state count. Pending: do the per-region leaves of one conformer MERGE
to one state (fixed 12 A vs auto floor; merge_rule_test on synth_null_50k / synth_two_50k)? 700k chain: family 1 done
(25 states, 275k), family 2 running; 4M merging its 410 GB stack since 22:13. Real data region_1 (clean): seeded 6
states, reproducibility 0.915 (random pending).

**Merge rule on the null / known-K sets (2026-10-19 00:xx):** null_50k: fixed 12 A floor -> 0 merges, 12 states of
1-7k (the ONE conformer sliced by the region boards; region leaves differ reproducibly at 12 A by subset systematics:
junk share, views) = over-split; auto floor -> 2 states (27k + 10k). two_50k: fixed -> 8 states, purity 0.97, both
conformers covered (4 slices per conformer); auto -> 2 states but purity 0.77/0.81 with ONE majority conformer =
over-merged: a region-A leaf (pure 18) and a region-B leaf (pure 49) compared inside A u B miss the arm when neither
region covers it. Fix: merge only if the difference is unresolved both inside the union mask AND over the whole
envelope (max of the two FSCs); rerunning on null_50k, two_50k, 55_200k (merge_rule_test). Target: null -> 1 state,
two -> 2 pure states.
null_200k (2026-10-19 01:10): every region board -> ONE supported leaf, purity 1.00 (11 boards): no over-splitting at
the support level at 200k either; the merge across regions is the only place a homogeneous set can end with > 1 state.
Band-scan auto rule (best 3-shell band from 40 A down to the pooled resolution, inside regions AND envelope) under test
on null_50k / two_50k / 55_200k / null_200k.
Band-scan rule result (01:4x): merges NOTHING on the null sets (15 / 14 pure-conf49 states, junk 0.00) and only 1 of
13 on 55_200k: same-conformer region leaves differ reproducibly in SOME low-resolution band although they hold the same
conformer and no junk -> the leaves are VIEW-SORTED (on homogeneous data the boards sort by view coverage, and
view-biased subsets give reproducible map differences in both halves; cf. view-imbalance-inflates-between-state-metric).
Tests: view_sort_check.py (direction histograms per leaf vs whole set) and VIEWMATCH=1 in consolidate_states.py (thin
every state to the whole set's direction histogram before reconstruction), merge tests rerun with it on null_200k and
two_50k. The merge rule must compare view-matched reconstructions.
View sorting confirmed (view_sort_check, 02:0x): leaf direction histograms vs the whole set: KL 0.02-0.12 on null_200k
(random subsets of the same size ~0.002), max bin ratio 1.4-2.3; two_50k KL 0.016 vs 0.005; 55_200k KL 0.05 (one tiny
leaf KL 1.4). Mild but 10-20x a random subset: enough to make two same-conformer leaves differ reproducibly at low
resolution. View-matched consolidation (VIEWMATCH=1: thin each state to the whole set's histogram) under test on
null_200k / two_50k (jobs 53995221/2).
TRAP (02:3x): VIEWMATCH in consolidate_states.py thins the MERGED states, i.e. after the merge; the merge itself reads
the chain's leaf half maps, so the first "view-matched" merge test changed nothing (0 merges again). Proper test =
vm_merge_test.sbatch: identity consolidation with VIEWMATCH=1 -> view-matched leaf half maps (consolidated_12_vm) ->
merge_data_derived.py on those (envelope mask) -> purity (jobs 53995274/5 on null_200k / two_50k).
