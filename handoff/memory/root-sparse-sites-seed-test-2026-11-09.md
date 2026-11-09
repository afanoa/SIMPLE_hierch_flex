---
name: root-sparse-sites-seed-test-2026-11-09
description: "VERDICT: the sparse pursuit (ladder + noise model, 7 seeds, same 50k subsample of the 692k root) reproduces the candidate REGIONS (top three = lower F1 / upper F2 / central F4 in 7/7) but makes one site everywhere and flips F1 vs F2 (4/3; gap < 2 SE in 4); fails the 6/7 leading-block bar; closed lives in F1 (= old site 1, never chosen): LDA .72-.77, label-free poles x1.5-1.8; f2 is no pursuit block"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T17:25:17.766Z
---

Coordinator task (user: "Yes, test it"), 2026-11-09.
- Notes: `hetero_tree_20261028/handoff/ole-panel/ROOT_SPARSE_SITES.md`.
- Runs: `ole_panel/rootsp2/s1..7`, `read.txt`.
- Cost: 67.6 core-h, of which 3.6 went to a failed first launch.

**Regions are reproducible.**
- The top three block families are the same in all seven seeds:
  - F1, lower, at (30, -43, -41). This is the old layer's site 1, 10 A away: the top-variance candidate the old site_rank never chose (margins 3.2-20.1).
  - F2, upper, at (3, -21, 42). This is the old region 9 (whole / f1), and the r131 closed block sits 9-10 A away.
  - F4, central.
- The site supports agree at Jaccard >= .59, but there is one site per seed, covering everything.

**The leading block flips.** F1 leads in 4 seeds, F2 in 3. The top-two gap is 0.7-1.5 SE in 4 seeds and 5.6-7.4 SE (F1) in 3. The pre-set bar (>= 6/7) fails at block level and is met only trivially at site level.

**Closed signal.**
- The root's closed signal is in F1 in every seed: the best block's LDA is .72-.77, and label-free poles reach x1.5-1.8. On 692k (pr692 block 3) it reads .802.
- F2 reads only .55-.62 at the root, against .68 at r131.
- The single site's label-free split gives x1.17-1.36, against x1.3 for the f2 starting groups.

**The old f2** appears only as late, weak blocks (gain rank 15-20) and has no counterpart in the pursuit.

**Why:** the old layer's seed dependence ([[factor-mask-split-is-state-set-region-2026-11-09]]) motivated replacing it.

**How to apply:**
- The pursuit fixes the candidate set, not the choice. Under the 2 SE rule, keep F1 and F2 as tied.
- F1 is where the root's closed signal is.
- The pursuit has no flex-fit input: sparse_sites skips the embedding, and SIMPLE_SEED drives both the preparation sampling and the pursuit.
- Launch traps: [[flex-pca-subset-is-pindfile-cap-needs-simple-path]].
