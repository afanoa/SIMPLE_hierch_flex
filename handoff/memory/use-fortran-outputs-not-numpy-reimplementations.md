---
name: use-fortran-outputs-not-numpy-reimplementations
description: "user 2026-10-13 ('lol why are you using a numpy prototype, the code is somewhere'): analyse real data with the existing Fortran code or its outputs, never a numpy re-implementation of the model"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T14:25:04.784Z
---

When a question is about what an existing SIMPLE program (cls_expansion, flex_pca) does on real data, use that
program's outputs (labels, averages, coordinates) or run the program; do not re-implement its model in numpy. The
numpy per-class PPCA on RNase P disagreed with the Fortran expansion the user had already looked at, and the user
caught it immediately ("this seems incorrect", "that's better" once I switched to the Fortran labels).

**Why:** the Fortran code carries the masking, solvent prior, bands and nuisance handling that the prototype lacks;
a mismatch then reads as a scientific result when it is an implementation gap.
**How to apply:** first look for the program output or a dump switch in the trees (and which build has it); if the
needed output only exists in another branch, say so and get that branch. Python toys stay for NEW model ideas with
known truth ([[prototype-in-python-before-fortran]]), not for re-deriving what Fortran already computes on real data.
