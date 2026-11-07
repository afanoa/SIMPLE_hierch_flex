---
name: site-mask-rework-direction-2026-11-07
description: "USER DIRECTION 2026-11-07: site ranking + masks are the pipeline's weakest part; prove the tree fails on a hard (IgM-regime) test, then explore flex prior, depth wiener, global coordinate there; formulation not cohesive, throws information away. Plan = guide section 7.7 (T1 focused per-site flex on the IgM fixture)"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T06:01:35.391Z
---

User, 2026-11-07, after the overnight verdicts: "I think we need a harder test (i.e. like the IgM regime, where can prove the
current tree fails) ... flex prior, depth wiener, and global coordinate needs to be explored. Of course there needs to be tangible
benefit for this stuff, but the current formulation is not really very cohesive and throws away a lot", then "perhaps we need to
address this whole site ranking and mask stuff. And perhaps the work we have done can have a use here. Because to me that seems to
be the weakest part of the pipeline". Same evening: do not push the manuscript drafts to the fork.

Plan of record: `doc/implementation_notes/hetero_tree_testing_guide_2026_11_07.md` section 7.7 on fork/hetero-tree. The hard test is
the IgM fixture (`igm_cryolo/tree/fixture_igm`: current tree = one leaf, a1 found only by islands + xval site stage B + GAPw at NMI
.146) with mm3q and qtoyflex as second tier and the section-3 gates as guards. First decisive test T1: focused flex_pca
(`pcg_mskfile` = each arm island site, `arm_masks.mrc` as control) on fix and null at the root's poses.

**Why:** every frozen-pose coupling failed on fixtures the tree already handles; the oracle put the headroom in seeds/references,
which the site/mask/ranking chain produces.
**How to apply:** judge flex prior, depth wiener and the global coordinate by benefit on the IgM fixture (null twin unsplit), not
on qtoy2/r131 alone; keep the r131 corner gate. Related: [[igm-fixture-xval-stageb-gapw-2026-11-06]],
[[igm-root-mask-is-flex-solvent-imprint-2026-11-05]], [[stageb-oracle-references-limit-2026-11-07]],
[[hetero-tree-algorithmic-audit-2026-11-01]], [[mac-retired-handoff-2026-11-07]].
