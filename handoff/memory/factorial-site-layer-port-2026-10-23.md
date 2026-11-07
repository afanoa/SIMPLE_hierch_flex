---
name: factorial-site-layer-port-2026-10-23
description: "lean factorial per-site layer built into hetero_tree (uncommitted, m8c6-m8c14); toy gate variants V0-V4, EPM labels, V3 chosen for full692; p27 toy validation; TRAP sacct COMPLETED != tree finished"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T03:53:16.899Z
---

Built 2026-10-22/23 night in worktree sharp-proskuriakova-b68280 (UNCOMMITTED; HEAD 17b2ea215), relayed by the coordinator
session; batch2 ~/SIMPLE_m8c builds 6-14, snapshots ~/SIMPLE_stree/snap/m8c6..m8c14 (each gated, fast gate exit 0).
- site_rank: nsample = fit on a random subset then score everyone; site_size=yes sizes given rankings (no EM); curve files carry
  hold-out amplitude + equivalent pure members (EPM); site_neff.txt + Fisher SE per endpoint.
- hetero_tree site_layer=yes (default): every candidate ranked; copies/dedupe by Spearman >= 0.5 only (no Dice); climbs decided on
  pole-1 Spearman alone; refinement (pole-set classes, falloff_mask) with guard = no hard-share flip + peak non-inferiority (1
  combined SE); sizing step; labels = top-EPM (majority pole = complement of minority poles, m8c14); maps = top-N*.
- TEST SWITCH env SIMPLE_FSITE_GATE=V0..V4 (outputs suffixed _vN) still in the code: hard-wire the chosen gate before committing.
Toy fulltoy_cdisc (factorial session scores): V0 misses closed (best AUC 0.627, 3 false sites); V1 0.816 / 4 false; V2 0.627;
V4 0.816 / 3 false; V3 (refine-first narrow fails) recovers via r1_s4 (0.789 at root) with 3 false. EPM fixes label inflation
(r_s2 closed pole 101k -> 10k labels, precision 0.445 vs Python 0.385). Coordinator rule -> full692 runs V3 + EPM (job 54110499,
started 2026-10-23 ~00:00, nodepar 6).
**Why:** the 3-SE halves gate drops the rare closed site (n_eff 50-110 masks); refine-first adds power without admitting nulls.
full692 V3 DONE 2026-10-23 04:36 (factorial_v3/, final labelling rule): 188 candidates -> 88 kept (39 home, 49 pooled),
31 copies, 33 dups, 36 no pole; 516 core-h / 4.6 h wall (~4x the 110-130 estimate: 145 pooling climbs at depth 0-1 at the 100k
fit cap = 199 core-h; candidates climb BEFORE dedup). First attempt (m8c13) crashed: a 1-particle rest class had no map ->
fixed in m8c18 (with_maps skips classes without maps).
**How to apply:** price per-site passes from climbs x (1.8 core-h at the root, 1.0 at depth 1), not per site. TRAPS: score_port_sites
--table reads only the first/last pole column (port sites can have 3-4 poles); the flip clause is decisive (a collapse raised a
null site's hold-out 0.28 -> 0.64). Discovery re-cut dry run: twin a r2 re-cut recovers 4 of twin b's 6 sites (Dice 0.46-0.77),
awaiting the user's OK. Related: [[state-tree-mergeup-port-2026-10-21]], [[port-tree-site-layer-recovers-closed-2026-10-22]].
p27 toy validation (m8c27, validate_tree.py, table in factorial_site_layer_results_2026_10_21.md): purity 0.72-0.73, domain A split
on every cdisc path, cdisc closed site seed-dependent (cd2 AUC 0.935, cd1 none). TRAP: a toytree job ends COMPLETED in sacct even when
hetero_tree failed (p27igm: End of file at simple_hetero_tree.f90:3203, the short read fixed in 387e6d08d) -> wait on hetero_tree_labels.txt.
