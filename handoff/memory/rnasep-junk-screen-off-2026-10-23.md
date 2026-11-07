---
name: rnasep-junk-screen-off-2026-10-23
description: "USER DECISION 2026-10-23 (via coordinator): junk screen OFF for RNase P - every hetero_tree / site_rank run on RNase P data (incl. the full-data tree) passes junk_screen=no; build m8c33 (HEAD b6db9a95f); send the launch line to the coordinator first"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-10-25T12:58:14.515Z
---

The user looked at the low-"amplitude" group that the junk screen flags at the full692 root: it looks like different RNA
conformations, not junk. Low brightness against the mean map also marks conformational outliers.

**Why:** the screen removes real conformers of a flexible RNA.

**How to apply:** any `prg=hetero_tree` or `site_rank` on RNase P data gets `junk_screen=no`. Launch from the current tree
build (m8c33, HEAD b6db9a95f at 2026-10-23: one pass per stage, node_rows=yes, oneflex, sta0_npeaks=2 by default) and send the
launch line to the coordinator session before submitting.
- FULL-DATA tree (4.5 M, user-approved 2026-10-23): NOT on m8c33. Wait for the port's snapshot with node-sized node
  projects + memory scaled for site/split (~1.5 days). Recipe ~nparts 150 nthr 8 nodepar 8, workers "-p csbdevel,norm",
  master 16 cores / 96 GB, junk_screen=no; target 8-10 h.
- The code default stays empty-only (the user said "for rnase p", not globally).
- hetero_tree forwards junk_screen to both site_rank calls (the node's S_SITE and the per-site fs_* layer), so no node gets a
  site_junk.txt or junk labels. A resumed tree could still read a stale one left in a seed dir.
- The contrast fit still runs, so stage A/B contrast stays on.
- TRAP (2026-10-25): junk_screen=no never reaches the split test's own screen. That is state_mergeup%screen, called unconditionally
  from test() in state_split and the merge-up; it flags a unit whose top-band per-particle SSNR vs the consensus is < JUNK_SSNR 0.35.
  - At the full692_design_s1 root it junked 124,999 particles (leaf r_s2), which were the molecule: corr 0.956 with the consensus,
    FSC 0.76 at 8 A to the kept unit, 58 % of the difference inside site 6's mask, and site 6's pole extreme.
  - Fix approved by the coordinator (2026-10-25), build m8c80t: forward junk_screen to the split and the merge-up, report-only
    under 'no'; the default stays yes.
  - Every RNase P tree before m8c80t ran with this split-step screen ON. Check a tree's build before trusting its junk leaves.
- Never edit a script a running job uses.
- 2026-10-23 15:30 (user): the full692 twins (m8c20, old rule) were CANCELLED by job id and moved aside to
  stree_port/full692_tw{1,2}_cancelled_m8c20. ONE replacement tree runs: stree_port/full692_m33_s1, job 54157123,
  stree_f692m33.sbatch (m8c33, SEED=1, junk_screen=no, site_layer=no). Estimate 550-800 core-h, 7-11 h.
- r3 was called "quality" from a lower refinement score, but conformational outliers also score lower, and 54 % of the group
  sits in r3 (maps: stree_port/rnp_root_junkgroup_maps).
Related: [[tree-contrast-and-identity-junk-2026-10-23]], [[junk-detected-per-node-not-only-root-2026-10-22]],
[[rnasep-full81k-preprocessing-run-2026-10-22]].
