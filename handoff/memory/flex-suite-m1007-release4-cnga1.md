---
name: flex-suite-m1007-release4-cnga1
description: "2026-10-07 suite tag m1007 = origin/master bd3997758 (hael's flex refactor pass 1): ribo/splice/trpm4 reproduce mix16 inside the run-to-run band (no regression); m1007 is the master baseline now; cnga1 blocked on release 4 (pruned stack, no indstk) and dropped by the user"
metadata:
  node_type: memory
  type: project
  originSessionId: 1a1d7175-9257-486d-a2bb-2f28c38b4e51
  modified: 2026-10-07T16:09:13.777Z
---

**Verdict:** latest master (bd3997758, includes hael's "flex AI slop refactoring, pass 1" 878adf2d4 and
"Release 4" d1b92e473) passes the flex_pca suite: ribo, splice and trpm4 match mix16 (2026-09-18,
unified p5 build) inside the run-to-run band, no errors. `m1007` is the reference tag for the next
suite run on master (`python3 ~/flex_verify/tools/suite_score.py unified6 mix16 m1007 <newtag>`, via
srun on csbdevel, never the login node). The user looked at the splice state maps: "looks fine".

| ds | tag | a | K | ge.7 | ge.85 | meanmin | states: occupancy / eo-FSC(0.143) | wall |
|---|---|---|---|---|---|---|---|---|
| ribo | mix16 | .43 | 4 | 4 | 4 | .560 | 8 @ 4.2-4.9 A | 28 min |
| ribo | m1007 | .43 | 4 | 4 | 4 | .557 | 10026/9684/14331/16292/15627/16019/14450/8818 / 4.8/4.6/4.3/4.3/4.2/4.2/4.3/4.9 | 27 min |
| splice | mix16 | .41 | 5 | 5 | 2 | .642 | 5 @ 8.9-9.7 A (unified6 meanmin .664) | 22 min |
| splice | m1007 | .41 | 6 | 6 | 2 | .627 | 13952/9941/6714/7168/9008 / 8.8/9.4/9.5/9.7/9.7 | 23 min |
| trpm4 | mix16 | .40 | 9 | 9 | 6 | .848 | 5 @ 8.0-8.5 A | 37 min |
| trpm4 | m1007 | .40 | 12 | 8 | 5 | .839 | 61070/31547/28980/25813/14916 / 8.0/8.1/8.3/8.1/8.3 | 42 min |

Reading it: splice now hits the scorer's reference row exactly (ge.7 6, K 6); its merged "min matched
|cos|" drop (.54/.42 -> .39/.30) is one tail component, meanmin moved less than unified6-vs-mix16.
trpm4's both-halves count is known noisy (4 to 9 across earlier runs); the only thing to watch is its
wall time, +13% (42 vs 37 min). The run log no longer prints the "sigma2 groups ... using the pooled
(global) noise" line on splice, but both runs used the same SIGMA2 GLOBAL BOOTSTRAP (message removed,
behaviour unchanged).

**Setup:** batch2 worktree `~/SIMPLE_m1007` (detached at bd3997758, login-node build, system OpenBLAS,
~10 min, make rc=0). Jobs 53718712 ribo / 53718713 splice / 53718714 trpm4; run dirs
`<dataset root>/flex_suite/m1007/<ds>` (ribo root `10028/10028_simple`, splice `10180`, trpm4
`trpm4_flex_afan`). suite.sh / suite_arm.sh now take MASTER_MEM / PART_MEM (MB); this run used
64000 / 32000 with no OOM (measured peaks 10-34 GB masters, 8-17 GB workers), and the arm exports the
gcc15 LD_LIBRARY_PATH. Local mirrors in scripts/flex_pca_verify (untracked); nothing committed.

**cnga1 is blocked on release 4:** `map_ptcl_ind2stk_ind` stops with "nptcls_stk not present in stack
field". `fix_projfile` (validate_projfile renamed; repairs IN PLACE, so run it on a copy) sets
nptcls_stk=274333 from the header but refuses indstk: 245,424 project rows for 274,333 images (pruned
without rewriting the stack), so the row-fromp+1 mapping every earlier suite run used is unprovable.
Failed run dir `cnga1/flex_suite/m1007/cnga1` (the arm refuses reuse; rename before a rerun);
unchanged copy `10_refine3D_auto/cnga1_test_r4.simple`. User's decision 2026-10-07: drop cnga1 for
now; the suite is ribo/splice/trpm4 until the project is re-extracted or the user picks the
legacy-mapping route.

**How to apply:** compare the next master suite run against m1007 (same build recipe, same memory
overrides). Any older dataset project may hit the release-4 wall; check with fix_projfile on a copy
(srun on csbdevel) before submitting.
Related: [[flex-pca-mcfa-ab-suite-2026-09-18]], [[cnga1-is-the-apo-closed-dataset]],
[[size-job-memory-to-measured-rss]].
