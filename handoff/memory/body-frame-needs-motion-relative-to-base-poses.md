---
name: body-frame-needs-motion-relative-to-base-poses
description: "2026-10-03 toyJ t3s0.02: body-frame composition on refine3D poses LOSES (true motion 0.95->0.65); steering domain coefficients lose too; toy3's centred hinge lets refine3D track the domain"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-03T21:14:44.526Z
---

Arms on batch2 `~/toyJ/toyj_bodytab.py` + `toyj_body.sbatch` (env TAG BASE MOT; MOT = none | gt | steer:<dir>:<scale>[:rank1]),
scores in `t3s0.02/body/SCORES.txt`. Region-1 body FSC@8/6 A vs the exact piece; the core is scored in reg96_2.

| base poses | none | true motion | steering x1 rank1 | steering x2 |
|---|---|---|---|---|
| true | 0.63/0.56 | 0.97/0.95 (core 0.73) | - | - |
| refine3D est | **0.95/0.90** (core 0.96/0.88) | 0.65/0.54 | - | 0.68/0.51 |
| robust frame r1 | 0.78/0.70 (core 0.96) | 0.90/0.79 | - | 0.61/0.38 |
| steered | 0.92/0.85 | - | 0.89/0.80 | 0.75/0.61 (rank1) |

- **Why:** refine3D poses already absorb most of the domain motion, so a motion composed on top of them must be the residual relative to those poses, not the absolute conformation.
- **Toy caveat:** toy3's hinge is at the box centre (−0.2, −0.5, −0.4), so a global rotation exactly mimics the domain motion and refine3D can track the domain. The fair test is an off-centre hinge (`toy3_offc`, HINGE_MODE=outer).
- **Steering coefficients:**
  - The raw coefficients carry 5–11° off-axis noise. A rank-1 SVD restriction recovers the hinge axis (|cos| 0.98).
  - Even so, per-particle corr 0.75 is too noisy to beat the consensus.

Related: [[toy3-body-frame-test-system]], [[robust-frame-loop-removes-heterogeneity-pose-bias]], [[joint-pose-heterogeneity-toy-verdict]].
