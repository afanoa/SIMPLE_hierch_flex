---
name: null-safe-joint-seeding-replicates-at-r1-2026-11-09
description: joint seeding (null-calibrated mixture over the kept sites' latent axes, K* 6) refined through the tree's own calls at r1 replicates across seeds (ARI .67-.69, 3/3 hits x2.3-2.5, six children) where single-site seeding is a coin flip (2/6, ARI .31 across basins); costs recall (.23 vs .38)
metadata:
  type: project
---

node_levers/nl_joint2.py + nl_jseed.sbatch (2026-11-09). Null: one Gaussian with the node covariance, each particle's deviation
scaled by its noise level from residual and contrast; stepwise K accepted above the null's 95th percentile. Toy nulls K* 1;
micnull K* 3 but the tree's locality guard folds its refinement to one leaf; qtoy2 K* 6 -> 6 of 8 conformers at purity >= .5
(tree 4 of 8), holo lost; qtoyflex K* 2 -> worse than the tree (1 of 3). r1 (arm a replays, seeds 1-3): K* 6, identical groups
across seeds (flex latent seed-independent); refinement keeps six children of 32-49k (sizes within 1 %), closed-leaning child
44.7k at .045-.048 (x2.34-2.47), reference recall .23-.24; cross-seed ARI .67-.69 vs .31 across single-site basins; production
r1 -> r13 was 80,775 at .041 (x2.12), recall .38. ~60 core-h for the three r1 runs.

**Why:** the user's rules: diversity first, generalize across datasets, not replicable = redesign. Joint seeding is the first r1
seeding that replicates; single-site seeding's hit / miss was set at the site layer ([[node-levers-replay-verdict-2026-11-08]]).

**How to apply:** candidate replacement for the site layer's single-pole seeding; before adopting, settle the recall cost (the
closed rows spread over three of the six children), the holo loss on the toys, and run it on the panel (IgM fix root, OLE rerun)
as the coordinator planned. Built with [[wip-settling-pass-prob-beta-cold]] unaffected (seeding, not refinement).
