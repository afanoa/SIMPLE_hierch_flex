---
name: trpm4-8extract-pixel-size-defect
description: "trpm4/8_extract has smpd mislabelled 1.3 when the true value is ~2.35 A/px (1.807x); use trpm4_correct instead, and sanity-check particle extent against the known molecule size before any campaign"
metadata:
  type: project
---

Found 2026-09-20 after a whole trpm4 campaign produced only "breathing" states.

**The defect.** `/mnt/beegfs/elmlund/processed_test_datasets/trpm4/8_extract` and
`.../trpm4_correct` process the SAME movies but motion-corrected to different
micrograph dimensions, and BOTH record `smpd=1.3`:

| | micrograph | smpd label | implied field |
|---|---|---|---|
| trpm4_correct (right) | 2306 x 2306 | 1.3 | 2998 A |
| trpm4/8_extract (wrong)| 1276 x 1276 | 1.3 | 1659 A |

Same detector, so the field cannot differ: 8_extract's true sampling is
**2998/1276 = 2.35 A/px**, i.e. the label is off by **1.807x**.

**Independent confirmation:** measured class-average extent of the TRPM4 particle is
**157-160 A in trpm4_correct vs 86-91 A in 8_extract** (ratio 1.82 == 2306/1276 =
1.807). TRPM4 is a ~150 A channel, so trpm4_correct is the correct one.

**Consequence:** everything built on 8_extract has resolutions inflated 1.807x
(reported 3.44 A really ~6.2 A) and every Angstrom-denominated knob meant something
else (`lp=14` was really ~25 A; that is why flex delivered breathing rather than
states). `kstop` was 16 in band vs 27 on the reference at the same nominal lp.

**Use `trpm4_correct` instead.** `19_abinitio3D/trpm4_correct.simple` = 162,326 ptcls,
box 300 @ 1.3, mskdiam 250, pgrp c4, poses already there; it is the project behind the
good 2026-09-16 five-state result and reproduces distinct states on current master
(user-confirmed 2026-09-20). `22_abinitio3D` is a 137,693-ptcl subset.

**The lesson:** measure the particle extent in the class averages against the known
molecular size BEFORE running a campaign. An 88 A extent for a ~150 A channel should
have been caught on day one; instead it cost a full sieve + 2D + two abinitio3D
branches + three flex arms. `print_project_field oritype=mic | grep xdim` next to the
smpd is a two-second check.

Related: [[trpm4-flipqr-sieve-campaign]], [[flex-pca-conventions-branch-state]].
