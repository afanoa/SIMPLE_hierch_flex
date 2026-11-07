---
name: cnga1-continuum-is-the-micelle-halo
description: "2026-09-14: cnga1's reproducible heterogeneity is the low-density peripheral halo (detergent/micelle), NOT the channel gate -- difference power is enriched 3.5-4.3x at r=50-80 A and DEPLETED 0.44-0.57x in the pore (r<30 A); mskdiam/pcg_mskfile is the lever"
metadata:
  type: project
---

**Measurement** (states_ax0, the 6 well-populated states, box 288 at 0.844 A/px, difference power
= sum over states of (state - ensemble mean)^2, normalised by the consensus power in the same shell):

| cylindrical radius from the C4 axis | % diff power | % consensus power | enrichment |
| 0-10 / 10-20 / 20-30 A | 2.5 / 8.8 / 17.9 | 4.5 / 20.2 / 31.7 | **0.56 / 0.44 / 0.57x** |
| 30-40 A | 18.9 | 22.7 | 0.84x |
| 40-50 A | 16.3 | 11.6 | 1.40x |
| 50-60 A | 21.2 | 6.1 | **3.50x** |
| 60-80 A | 14.3 | 3.3 | **4.31x** |

35.5 % of the variability sits beyond r=50 A where only 9.4 % of the density lives, and the pore
region is *depleted*. Axially it is mildly mid-plane concentrated (1.36x at |z|<10 A).

**Geometry**: the consensus peaks at r=10-30 A (0.544/0.548 mean |rho|) and decays to 0.146 / 0.081
/ 0.024 at 50-60 / 60-70 / 70-80 A. `automask3D_state01.mrc` -- the `pcg_mskfile` every cnga1 arm
used -- reaches **r = 77.8 A**, so the estimator's support encloses the whole micelle. The protein
is r <~ 40 A. This is the same failure mode as [[flex-pca-solvent-mask-dominates-rna-basis]] (there
mskdiam=210 on a ~140 A molecule), and the lever is the same: a protein-only support.

**Do NOT test this with the automask**: it is 7.98 % of the box and encloses the belt, so
"99.3 % of the difference power is inside the mask" is vacuous. Separate radially.

**Why everything else was consistent with this**: a disordered belt has many near-degenerate
deformation modes, so the spectrum stays flat (max/min 1.35-1.59) and the reproducible rank keeps
growing with the budget (8 of 10 -> 10 of 16 -> 15 of 24, mean |cos| 0.798/0.738/0.737) without
saturating, and a micelle is a continuum by construction so the dip statistic is 0.
Related: [[cnga1-is-a-continuum-not-discrete-states]], [[flex-pca-states-stratify-by-neff]].

---

**MASK PROVENANCE + THE envmsklp LEVER IS REFUTED (2026-09-14, measured).**

`pcg_mskfile=10_refine3D_auto/automask3D_state01.mrc` is NOT a solvent-support mask despite the
name: it is the **envelope-FSC** mask, written on the `l_envfsc` path
(`simple_rec3D_pcg_strategy.f90:280`) by `automask3D(..., lp_override=params%envmsklp)`
(`simple_halfmap_diagnostics.f90:35`), i.e. the consensus is low-passed (`self%bp(0.,amsklp)`)
before thresholding. The cnga1 run recorded **`envmsklp=20.0`** = `ENVMSKLP_DEFAULT`
(`simple_defs.f90:90`) -- **not `amsklp` (8.0), which `lp_override` replaces.** The sibling file
`nu_envmask3D_state01.mrc` is the different NU *evidence* envelope (`nu_msk_sig=3.0` MADs,
`amsklp=8`, +1 A grow, 6 A edge, 0.1 min-vol CC); refine3D's own PCG deliberately does NOT use it,
so passing automask3D was right by convention.

**Regenerating at amsklp=8 vs 20 changes almost nothing, and the "20 A pulls the belt in" mechanism
is WRONG.** `simple_exec prg=automask vol1=rec_final_state01_pproc.mrc smpd=0.844 mskdiam=180
amsklp=<lp> nthr=8 mkdir=no`, run on csbdevel (the opt8 binary SIGILLs on the batch2 transfer node
-- `-march` mismatch; always sbatch to csbdevel). Control reproduces the archived mask at
**IoU 0.921**.

| shell | archived (lp20) | new lp20 | new lp8 |
| 0-20 A | 1.000 | 1.000 | 0.969 |
| 30-40 | 0.995 | 0.998 | 0.957 |
| 40-50 | 0.971 | 0.975 | 0.878 |
| 50-60 | 0.800 | 0.799 | 0.720 |
| 60-70 | 0.476 | 0.435 | 0.351 |
| 70-80 | 0.102 | 0.057 | 0.049 |
| frac / maxR / medR | .0798 / 79.1 / 52.0 | .0761 / 77.0 / 51.0 | .0682 / 77.8 / 50.2 |

lp8 removes only **10.4 %** of mask volume and **23.3 %** of the archived mask volume in 50-80 A,
leaves 72 % occupancy at 50-60 A, and **max radius does not shrink (77.0 -> 77.8 A)**. It also
ERODES the protein interior (1.000 -> 0.969 at 0-20 A) because a finer low-pass makes the density
granular and thresholding then drops interior voxels. **So the envelope's generosity is set by the
THRESHOLD, not the low-pass: the micelle really does carry supra-threshold density to ~70 A.**
A protein-only support needs an explicit hard radial cut (r <~ 45 A), not an `envmsklp` change.

**Unit trap I hit:** cnga1 box 288 is **smpd = 0.844 A** (cella 243.072), NOT 1.06277. Using the
wrong value overstated every mask radius by 26 % (I reported max 99.6 A and an 80-100 A shell at
0.220; truth is 79.1 A and 0.000). `cl_rad.py` used 0.844 correctly, so the r=50-80 A enrichment
result stands -- only the mask description was wrong. Always read smpd from cella/box.
