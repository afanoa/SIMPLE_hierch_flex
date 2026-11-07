---
name: igm-picking-audit-2026-11-05
description: IgM (EMPIAR-11077) picking is adequate - our template picker finds ~85-89 % of the IgM molecules any of 5 pickers found (off-centre duplicates removed); face-on IgM is rare for every picker (3.7 %); ~55 % of picks are junk
metadata:
  node_type: memory
  type: project
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-06T04:35:29.783Z
---

Picking audit, 2026-11-05, run in `igm_cryolo/audit/`. It used DJ's 1,000-micrograph test batch (at least 10 sessions).
- Five pickers were merged: ours (template picks with pickrefs_v3), a 3D-projection template picker, segdiam (diameter >= 150 A), a Gaussian blob picker, and DJ's picks.
- The merged set was extracted, run through 2D with 150 classes, and the classes were scored against our de novo maps (`pa_score.py`; IgM cutoff .6, calibrated on the user's core/keep/drop classes).

**Result**
- IgM classes held 48,040 of the 130,852 merged particles.
- Raw recall: ours 51 %, 3D-template 49 %, DJ 34 % (650 micrographs only), segdiam 32 %, Gaussian 19 %.
- Most "missed" IgM are our own particles picked off-centre. A side view is a rod about 350 A long, and the merge radius was only 50 A.
  - 59 % of misses have one of our picks within 150 A.
  - Only 5,492 (11.4 % of all IgM) have none of our picks within 250 A.
  - Corrected recall is therefore about 89 % (about 85 % with a 200 A cut).
- Face-on IgM is 3.7 % of all IgM for every picker, so the shortage of flat views comes from the sample, not from picking bias. Our picker still clearly misses about 27 % of face-on IgM, against about 11 % overall.
- About 55 % of our picks fall in non-IgM classes (ice, aggregates, discs). Classification removes these.

**Why:** the user asked whether picking is adequate and whether to repick with flat views.

**How to apply:** don't repick the full data. A union with the 3D-template picks would add only about 10 %, so that's optional. When estimating recall, merge picks of the same molecule at a radius larger than the particle's long axis, or the off-centre duplicates inflate the miss rate.

Related: [[igm-11077-consensus-poses-tree-ready-2026-11-04]], [[igm-retain-side-views]].
