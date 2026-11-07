---
name: stageb-oracle-references-limit-2026-11-07
description: "VERDICT 2026-11-07 (oracle refine=state_eval on qtoy2 root/r1/r4 + qtoyflex root): stage B is limited by its references (count at the root, purity in a close-pair split, no junk/minority reference), not by poses; noise ceiling reached only where references match the classes; the 32 % 'confidently wrong' are K=5 merges"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-07T05:04:07.009Z
---

Side session for the hetero_tree paper. Write-up and copies: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/stageb-oracle/VERDICT.md`.
Runs: `math_reworks/latent/oracle/<key>`. Scripts: `latent/orc_*`.
The arms are the node's stage-B call at refine=state_eval maxits 1. The poses are B or T; the maps are stage B's, coarse
truth (K_B), or one per truth class.

- **Poses are not the limit.**
  - qtoy2 stage-B poses: 2.2-2.6 deg median.
  - True poses with stage-B maps: root group accuracy 0.951 -> 0.947; r1 0.756 -> 0.758. They move 11-18 % of labels but
    fix as many as they break.
  - Poses add +1.5 to +4 points only once the references are pure.
  - qtoyflex's heavy pose tail is the exception: holo F1 +0.065 from poses alone.
- **The references are the limit.**
  - Root: K = 5 for 8 classes, so pairs merge. One reference per class at the same poses gives class accuracy
    0.674 -> 0.919 (ceiling 0.936). At its own K the root is 1.2 points from the ceiling.
  - r1 (1 vs 2): 80/20 references, a state_eval fixed point. Pure references: 0.756 -> 0.839 (ceiling 0.891).
  - r4: at the ceiling. Junk holds closed F1 at 0.896 (0.939 with a junk reference).
  - qtoyflex: no closed state, closed F1 0.17. Pure arm/closed/holo references at K = 3 give 0.455; with true poses 0.570.
- **The "32 % confidently wrong"** at the qtoy2 root is 27.8 % merged plus 4.7 % misassigned. 96.7 % of the >= .99 wrong are
  merged, with poses and snrq identical to the right ones. The true misassignments are unconfident, low snrq and pose-tailed.

**Why:** answers "where is the headroom" after every frozen-pose coupling failed. The headroom is in the references (the split
fixed point, junk and minority references), not in posteriors or poses.
**How to apply:**
- Push the paper's improvements at the split and the reference set. Do not prioritise pose refinement.
- Tricks:
  - Write true poses or labels into a .simple by byte-patching the ptcl3D slots (header 20 x 5 int64; segment 6, 64 float32
    per record; slots e1 8, e2 9, e3 10, state 23, x 27, y 30: `orc_patch.py`). No import, nothing else changes.
  - Rebuild maps through refine3D's own PCG/focus/ML pipeline from fixed labels: run state_eval with SIMPLE_STATE_PRIOR set
    to 0 on the label and -1e4 elsewhere, maxits 2.
  - Child-node frames rotate 1.1-1.3 deg from the truth: fit and remove the global offset before using true poses.

Related: [[flex-to-refinement-coupling-verdicts-2026-11-06]], [[mest-formulations-verdict-2026-11-06]],
[[flex-refinement-integration-verdicts-2026-11-06]].
