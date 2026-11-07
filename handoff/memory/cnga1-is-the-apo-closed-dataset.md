---
name: cnga1-is-the-apo-closed-dataset
description: "2026-09-14 DECISIVE: the cnga1 dataset IS the apo/CLOSED sample (7LFT / EMD-23306) from Xue et al. Neuron 2021, so open-vs-closed gating CANNOT be recovered from it -- only one state is present; this explains the continuum, the flat spectrum and the depleted pore variability"
metadata:
  type: project
---

**The user's actual goal was to reproduce the OPEN and CLOSED states. It is not achievable from this
dataset, for a data reason, not a method reason.**

**Ground truth** (Xue et al., *Neuron* 2021, 109:1302-1313.e4, "Structural mechanisms of gating and
selectivity of human rod CNGA1 channel"):
| PDB | EMDB | state | condition | res |
| 7LFT | EMD-23306 | **apo / CLOSED** | K+/Ca2+, no cGMP | 2.6 A |
| 7LFX | EMD-23308 | cGMP / OPEN | Na+/Ca2+ + cGMP | 3.1 A |
Maps at `/mnt/beegfs/elmlund/afan/cnga1/refs/emd_233{06,08}.map`.

**Identification, two independent lines:**
1. **Grid fingerprint.** Ours: box 288, cella 243.072, smpd **0.8440** -- an EXACT match to EMD-23306.
   EMD-23308 is cella 231.840, smpd 0.8050. Same processing pipeline as the apo entry.
2. **Rotation-invariant profile correlation:** ours vs apo **+0.862**, ours vs open +0.773, against an
   apo-vs-open separation of +0.851 (i.e. the 0.089 gap is meaningful on that scale).

**Why every method knob failed, now explained.** The population is one conformational state, so the
latent is a continuum (dip 0.00000 at d=10 AND d=16), the spectrum is flat, and the pore is the
LOWEST-variability shell (enrichment 0.51x). Nothing is gating.

**The disjointness measurement (the method-vs-data discriminator).** The open-closed difference is a
pore-axis event: **41.8 % of its power at r_cyl < 10 A, 70.5 % inside 20 A, only 3 % beyond 45 A.**
The maximum projection of any of the 10 fitted eigenvolumes onto that difference is
**|cos| = 0.157** (8 of 10 below 0.08). So the dataset's heterogeneity and the gating transition
occupy **disjoint regions of the molecule** -- our modes are peripheral, gating is axial.

**RETIRED as pointless for this goal:** symmetry expansion (gating is concerted and C4-symmetric, so
expansion buys nothing -- see [[cnga1-c4-pose-folding-scrambles-protomers]]); a focused small-mskdiam
run (would focus on a region with no signal); further box/band/rank/iteration work as a route to
discrete STATES.

**What would work:** (1) obtain the cGMP dataset and treat open-vs-closed as a TWO-dataset /
supervised problem where the ligand condition is the label; (2) a sub-saturating cGMP sample is the
only way it becomes genuine within-dataset heterogeneity; (3) inject the now-known gating difference
into these real images at known amplitude to calibrate flex_pca's detection floor for this specific
motion (machinery: [[pfcrt-injection-floor]]).

**METHOD TRAP, reusable: never compare C4 (or any Cn) maps by direct voxel correlation.** Independent
refinements settle at an ARBITRARY azimuthal angle about the symmetry axis, so direct cc gave 0.056
(ours vs apo) and even 0.099 apo-vs-open -- both meaningless. Compare **cylindrically averaged
rho(r_cyl, z) profiles**, which are invariant to rotation about the axis, searching z-flip and small
z-shifts. That took apo-vs-open from 0.099 to 0.851. A near-zero cc between two structures of the
same protein is a broken comparison, not a result.
Related: [[cnga1-is-a-continuum-not-discrete-states]], [[cnga1-continuum-is-the-micelle-halo]].
