---
name: tree-launch-fast-defaults-2026-10-25
description: "STANDING (user via coordinator, 2026-10-25): every hetero_tree I launch uses the fast settings: nodepar matched to the cores (16+ for big trees), workers -p csbdevel,norm WITH SBATCH_CONSTRAINT='x6342|x6740P', a build >= m8c90t (concurrent site ranking, no flex sort)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T23:08:11.577Z
---

2026-10-25, the user, after the speed work: "with speed yes I agree. But I have yet to see convincing evidence that we have beat the tree from yesterday."

- Every tree I launch from now on uses:
  - nodepar matched to the core budget: 16 or more for big trees (it only schedules, results are unchanged);
  - step jobs on csbdevel,norm, via QPART / SIMPLE_QSYS_PARTITION in the ENV prefix, never in --export (it cuts at commas);
    EXCEPTION: csbdevel only while the user's production jobs hold norm (2026-10-25: refine3D 54295680, 48 x 8 workers);
  - SBATCH_CONSTRAINT='x6342|x6740P' exported in the launcher. norm's x6150 Skylake nodes SIGILL our -march=native builds
    ([[batch2-archopt-sigill-login-node]]);
  - a build at m8c90t or later (concurrent site ranking, chosen-first motion, flex sort REVERTED). m8c75t-89t sort the flex
    components, which lost the RNase P closed state ([[flex-component-sort-moves-site-seeding]]); don't launch on them.
- Launchers that do this: stree_port/stree_f692t.sbatch (692k) and stree_wtoy4.sbatch (toys).
- The design stays frozen unless the quality toy ("qtoy", built by the factorial session) shows a benefit. The qtoy runs three trees:
  - (a) yesterday's tree: full692_m33_s1 ran on m8c33 for the root, m8c36 for the body and m8c61t for the end phase (merge-up,
    plain maps), capped at maxdepth 7;
  - (b) m8c87t with the design keys off;
  - (c) m8c87t with the design keys on.
  - (b), (c) carry the sort; (b') / (c') on m8c90t are the clean arms (added 2026-10-25 about 19:00).

**Why:** the user needs a fast test loop to see whether development helps at all.
**How to apply:** use these settings without being asked. Judge design changes on the qtoy trees, not on K* or margins.
Related: [[tree-runtime-goal-5-7h-at-700k]], [[tree-design-feature-freeze-2026-10-25]], [[m33-tree-wall-profile-2026-10-25]].
