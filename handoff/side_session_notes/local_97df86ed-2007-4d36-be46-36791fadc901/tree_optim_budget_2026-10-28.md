# hetero_tree time, CPU and disk levers (session tree-optim-budget, 2026-10-28, final)

Everything lives on beegfs under /mnt/beegfs/elmlund/afan/hetero_tree_20261028/; cluster_scratch was wiped. Nothing committed.
- **recipe12b/**: the production engineering build.
  - recipe9 + a1_disk + a2_sitebatch + a3_sitethreads + a4_plain (simple_hetero_tree.f90), + equiv_port_default_on (state_mergeup),
    + a5_strongest (strongest_pole default yes).
  - build_r12b, fast gate 15/15.
  - recipe13/14 (factorial session, route A) are built on it.
- **recipe15/**: recipe14 + recipe15_rowlocal_sitelayer.diff (row-local node projects with site_layer=yes; node_proj=yes by default).
  - build_r15, fast gate 15/15.
- **recipe12/**: gates, subtrees, the copy job, and analysis/ (node tables, sacct).
- **from_scratch/**: copies of recipe1-12 patches and scripts, recipe11/12 sources, the build_r9/r10/r11 installs and the passed
  recipe11 gates (9.0 GB, verified).
- Superseded: recipe12 (split threads changed), recipe12c (fallback, unused), recipe10/11 (early).

## Lever classes (the user's rule: "we cannot sacrifice science for speed")
- **Same computation, verified bitwise on the RNase P r221 subtree** (7 nodes, all stage-B maps byte-identical) and on the
  r1211 / qtoy null / OLE null gates:
  - disk cleanup;
  - rank all sites at once on nodes over 200k rows (inert under strongest_pole);
  - site-rank threads by node size;
  - job CPUs = threads;
  - plain maps as nodes finish (base only; off with site_layer=yes).
- **Same computation within numerical noise:** recipe15's row-local site layer (r111 route A md1: labels identical, decisions identical,
  flex at 1e-6).
- **Not equivalent, dropped:**
  - shared-memory small nodes. Changing the part decomposition changed site discovery (7 vs 4 sites at r3121); labels at ARI
    0.01-0.64. A reproducibility sensitivity in discovery, not only a speed result.
  - split/merge-up threads by size: 1 row moved at r221 flipped r2211's strongest pole (margin 6.40 -> 1.67).
  - 16k rows per part (untested decomposition change).
- Bug in recipe9: prune_partials deletes nothing. simple_list_files_regexp returns dir/name and the prune deletes dir//'/'//name.

## Gates on recipe12b (maxdepth 4; baselines = recipe9 + strongest_pole)
| gate | recipe12b | baseline |
|---|---|---|
| r1211 s1 | closed F1 0.616 | 0.616, ARI 1.000 |
| r1211 s2 | 0.555 | 0.555, ARI 1.000 |
| qtoy2 md4 | NMI 0.647; c101 0.82, c18 0.91, c28 0.95, c43 0.46, c49 0.60, c1 0.77 | md7 cut at 5: 0.626, the same, c1 0.65 |
| qtoynull | 0.669, labels identical | 0.669 (known fail) |
| olenull | 1 leaf, labels identical | 1 leaf |
| closed_r112_base | F1 0.100, 0 corner leaves | build-7 0.105, 0 |
| closed_r111_routeA (base mode) | F1 0.137, 0 corner leaves (0.94/0.55 nearest) | route A reaches the corner |
Cost at the same settings with identical results: -15 to -27 % core-h. Site ranking -65 %. Disk -82 % (subtree peak 13.2 -> 2.4 GB).

## Projections (not runs)
- Base, maxdepth 4, 692k (recipe9 levels 120/48/61/66/45 min): recipe9 ~6.3 h, ~650 core-h; recipe12b ~5.7 h at nodepar 8, ~5.2 h at 12.
- Route A, maxdepth 4, 692k (build-7 sst7 shape: 104 nodes, 1/6/13/28/56; spans 153/48/36/30/26 min):
  - recipe16 node work 7.4 h at nodepar 8, 5.85 h at 12, 5.2 h at 16;
  - ~890 core-h (build 7: 978);
  - disk ~140 GB with recipe15 (build 7 at md4 ~970 GB);
  - plus the factorial end phase (unmeasured at 692k; 1.4 h on OLE).

## Late 2026-10-28: recipe15b / recipe15c (on recipe15; for the factorial session's recipe17b lineage)
- recipe15b_reuse_prune.diff (recipe15b/, build_r15b, fast gate 15/15):
  - The site layer reuses the node path's strongest-pole batch ranks, when strongest_pole=yes, junk_screen=no and node rows
    <= SITE_NFIT (100k). Below 100k the two ranks are the same computation.
  - prune_node also runs with site_layer=yes. It keeps n.simple and the site, nsite and sites folders.
- recipe15c_states_delivered.diff (recipe15c/, build_r15c, fast gate 15/15):
  - final() writes 1_hetero_tree/STATES_DELIVERED (ISO timestamp) and a log line once the merge-up and final_maps are done.
  - The factorial phase keeps running; the exit is unchanged. No marker if the merge-up failed.
- Both diffs apply on recipe17b/source with offsets only (recipe17b = recipe17 + recipe15). Stacked dry-run: recipe15c/dryrun_r17b_2254.
  Sent to the factorial session as recipe17c (15b, md5 1ae3add2) and recipe17d (15c, md5 ca293bf0, cut on 15b, applied after 17c).
- Equivalence runs: r111 md1 recipe15/equiv/r111_md1_r15b and recipe15c/equiv/r111_md1_r15c, vs recipe15/equiv/r111_md1_r15.
  Compare with recipe15c/cmp_trees.sh; cost with recipe15c/run_cost.sh. cmp alone fails on absolute run paths: normalize them.
- Verdicts (both same computation; every step exit 0):
  - Byte-identical: labels, merged and sites tables, all 98 shared maps, seeds.txt.
  - Identical once paths are normalized: leaves, final_maps/choice, the factorial table and sites, leaf_masks.
  - 15b drops one state_seeds line: the seeding site against its own re-rank (rho 1.000, never added).
  - flex_pca text differs at 1e-9 between any two runs. One site_halves statistic varies between two runs of the same code; its
    decisions do not.
  - 15b cost: core-h 24.11 -> 22.35 (fsite 3.83 -> 1.71); disk peak 15.9 -> 10.2 GB, final 10.0 -> 5.8 GB; wall 76 -> 63 min.
  - 15c: STATES_DELIVERED at 58 min, exit at 63 min, log line present.
- The 13 min wall gain is a nodepar=1 effect. Every fs job counts against one global nodepar cap (fs_nrun), and fs_home_ranked
  holds seeding until the node's home ranks finish. So 5-9 home ranks ran one after another.
- Measured inputs for the 692k end phase:
  - On build-7's route A tree, nodes <= 100k finished their site-layer ranks before their own node-path ranks (median 4-8 min
    earlier). Build 7 ran at most 8 fs jobs at once, with 4-5 home ranks per node.
  - So at nodepar 12-16, reuse saves mainly CPU (52 of the 101 site-layer core-h at md4) and at most ~2 min of wall per level.
  - 692k merge-ups took 38-90 min at 32 threads (m33 90, base80 54, base90 38, base92 53; build-7 base killed at 79).
    The lineage now uses 48 threads.
  - Before the merge-up job: landscape 5 min, then 11 min of controller work (finish, final_project) on build-7 base.
  - The factorial runs concurrently with landscape and merge-up. So the marker saves only the factorial's excess over the merge-up.
  - Factorial site ranks take min(32, rows/20k) threads, up to nodepar at a time. At nodepar 16, root-level climbs could
    pass 400 CPUs during the end phase.

## Route A 692k projection, maxdepth 4, recipe17 line with 17c/17d (build-7 sst7 shape, 104 nodes; not a run)
| nodepar | node work | to states delivered | to factorial done | node CPUs peak / mean |
|---|---|---|---|---|
| 12 | 5.85 h | ~7.0 h (6.75-7.6) | 7.35-8.35 h | 229 / 138 |
| 16 | 5.20 h | ~6.4 h (6.1-6.95) | 6.7-7.7 h | 273 / 155 |
- After the node work, states are delivered in ~1.2 h (0.9-1.75): landscape, controller work, then the merge-up.
- The factorial phase takes 1.5-2.5 h after the node work (rough).
- Node work ~807 core-h with reuse (859 without).
- Disk peak ~90 GB if the measured -36 % carries over from the ~140 GB recipe15 estimate.
- Risk at nodepar 16: in the end phase, root-level factorial ranks (32 threads, up to 16 at once) plus the 48-thread merge-up can
  pass 400 CPUs.

## 2026-10-29: junk leaf for junk split groups (junkleaf/, on recipe17d)
- Why the IgM junk group became node r4: state_split's junk screen did not flag it.
  - Root unit 2 (49k particles, truth 73 % junk) had high-band rel. SSNR 0.64; the cut is JUNK_SSNR 0.35. It was reported only as
    "partial ... tested as composition".
  - Flagged units are already set aside as junk leaves (split_junk.txt -> add_leaf .true.).
  - Below the root, a pure-junk node cannot flag itself: its consensus is junk.
- No cut on the unit metrics separates it from real units: real qtoy2 and OLE units sit at the same relative SSNR.
- The signal that does separate it: the seeding site's junk screen (junk_screen=yes) writes each particle's amplitude-rule junk
  posterior (site_junk.txt column 6, not applied).
  - The seeding site comes from sites/chosen_site.txt, else seed_site.txt.
  - On three IgM trees, every split group with a mean posterior > 0.5 is >= 72 % junk: root g4 0.57, r1 g4 0.60.
  - Every genuine group is <= 0.29.
- recipe_junkleaf.diff (junkleaf/, build_jl, fast gate 15/15), 6 hunks in simple_hetero_tree.f90, md5 of the patched file 595fa84d:
  - In a passing split with junk_screen=yes, a group whose mean posterior is > 0.5 becomes junk leaf <node>_j<group> instead of a
    node.
  - add_leaf takes an optional particle selection and name.
- Inert under junk_screen=no, so qtoy2, both nulls, OLE and RNase P are unchanged by construction: their runs have no site_junk.txt.
- Verified: the IgM maxdepth-1 replay (junkleaf/igm_md1_jl) resumes a copy of the recipe17b root at its split. It logs "mean junk
  posterior 0.569 ... -> junk leaf r_j4" (27,642 particles, truth junk 1.00); r1-r3 stay nodes.
  - End to end: labels, merge-up (0 junk passed through; r_j4 stays its own unit), final maps and STATES_DELIVERED complete.
  - Exit 0, all step exits 0, 30 min, 35.6 core-h.
  - The diff's md5 is c5ef1fb4. It went to the factorial session as the next recipe17 increment.
- Saving on the recipe17b IgM tree:
  - 7 of 53 nodes (r4, r41, r42, r421, r422, r14, r141), 19.7 core-h of 175 so far (that tree is still running);
  - its junk is no longer delivered as 7 ordinary leaves.
  - Cost: r1 g4 holds 1,332 genuine particles (28 %), and they go to junk with it.
- Verdict (coordinator, 2026-10-29): the factorial session ran the full IgM tree on recipe17f. The rule fired once, at the root, and
  lost 69 genuine particles; it saved 11 % of core-h. Truth was slightly down, most likely from the r1 seeding near-tie. The rule stays
  out of the default, and the build is kept as is.

## Open
- The root's projects keep the input's 932k-row file size (write_cut copies, then rewrites segments in place).
- SIMPLE's slurm template writes --mail-type=FAIL into every step script, against the no-email rule.
