---
name: flex-deformation-model-falsified-on-igg
description: "The deformation-field (warp) alternative to flex_pca's linear basis was tested on IgG on 2026-08-05 and failed; a rotated-gradient control scored as well or better at every rank"
metadata: 
  node_type: memory
  type: project
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T20:17:30.428Z
---

On 2026-08-05 the Phase 0 test in
`doc/for_developers/ideas/flex_deformation_field_model.md` §9 was run on the
`/Users/atifao/igg_flexpca_probe` run (box_crop=64, smpd 6.0, lp=15, 48 comps,
5 states). **Result: negative.** The differences between flex_pca's IgG state
maps are not deformations of the consensus map.

Decisive numbers, states 4-5 (the most reproducible pair, ceiling
corr(d_even,d_odd) = +0.947): fitting a warp basis on the even halfset and
scoring the prediction against odd gives +0.161 at K=9 and +0.435 at K=33 —
while a **90-degree-rotated** `Vbar` gradient basis (same rank, wrong alignment)
gives +0.161 and +0.463. Separation is 0.00 to -0.06 at every rank from 3 to 408.

Both alternative explanations are closed:
- **Instrument works.** A synthetic 4.4 A-mean hinge at the same band, mask and
  ceiling separates +0.76 at K=9 (warp +0.925 vs rotated +0.165).
- **Not a large-motion failure.** Separation at K=9 is still +0.48 at a 22 A-mean
  / 40 A-max deformation. Linearization is far more forgiving than the lp/2
  argument suggests.

**FINAL STATE: the deformation model PASSES on IgG with the right basis.** With
elastic-network-model normal modes (120 density-placed beads, ANM Hessian,
Rc=4xNN=39 A, modes 7..66), measured on two FULLY DISJOINT 10k-particle
replicates of the endpoint contrast: separation **+0.144 +/- 0.001** (half-range;
bootstrap sd 0.023), roughly 2x any generic dictionary at equal or larger K
(scalar K=136: +0.080; vector3 K=90: +0.072). It holds at +0.144 across replicates
whose ceilings differ threefold (+0.337 vs +0.103), which is the best evidence it
tracks structure not noise. Scalar-momentum (LDDMM form, smoothing kernel
essential) beats the 3-vector form per parameter but is noisier than ENM.
Instrument: `scripts/warp_phase0_v4.py` (v3 superseded).

**SPECIFICITY CONFIRMED on Tomotwin-100.** Same protocol, GT species bins (20
proteins vs 20 different proteins): ceilings **+0.99**, and enm(60) separation
**-0.036 +/- 0.030**. Every basis lies in [-0.047, +0.031]. So the elastic basis
fires on IgG's WEAK conformational signal (+0.144, ceiling 0.10-0.34) and is silent
on a compositional difference 3-10x LARGER. This rules out the main alternative
explanation -- that 60 modes were just absorbing any reproducible structure.
Instrument: `scripts/warp_tw_enm.py`.

This null is CORRECT, not a weakness: a warp moves density and cannot create or
destroy it. DynaMight states the same limit for itself -- it "cannot model
compositional heterogeneity" because it refines no per-Gaussian occupancy, and is
"by design constrained to only model mass-conserving heterogeneity"; compositional
heterogeneity must be removed by discrete classification BEFORE running it
(Nat Methods 2024). If asked, note this plainly -- it is easy to misremember as
DynaMight handling compositional heterogeneity. The 2025 deformable-3D-Gaussian
follow-up (arXiv 2512.21599) adds the occupancy parameter and does handle both.
Corollary: a deformation flex path in SIMPLE must keep `flexmodel=linear` available
for compositional data. The local relion clone has ONLY the DynaMight launcher
(`scripts/python_dynamight.in`) and job-option strings; the model is in the
separately-installed `il` conda package, not readable from the repo.

**Critical fact about the IgG data: it has only ONE conformational contrast.**
Deflated ceilings over 5 trajectory bins show only b00-involving pairs carry
reproducible shape change; frames 20-99 are one indistinguishable population and
0-19 the other. So "run more bin pairs" does NOT buy replicates -- 6 of 10 pairs
have no signal to explain, and averaging over them is meaningless. Replicate by
disjoint particle subsets instead. This is also why b40 ~ b80 looked like a
fromp/top bug in an earlier pass; it is the data.

**Two earlier claims retracted:** the "+0.194" figure used 20k bins and is ~25 %
optimistic; and "ENM is the only basis that never flips sign" was read off two
pairs and is noise (over 10 pairs enm(6) is 7/10 positive, no better than generic).

Two traps in the ENM: assert **zero_modes == 6** (300 beads in a 3960-voxel mask
with Rc=1.9xNN fragments the network -> 12 zero modes, and a whole run was
discarded to this); and do NOT tune on the synthetic hinge -- vector3 wins there
(+0.757 vs enm +0.461) and loses on real data. Idealized-hinge performance
anti-predicts real-data performance.

**REVERSED earlier the same day.** Grouping the SAME particles by ground-truth MD
frame instead of by flex_pca's embedding (20k-particle bins, frames 000-019 vs
080-099, Vbar direction deflated out of d) gives separation **+0.117 / +0.122 /
+0.123** at K=9/33/90 — the warp basis explains 40 % of the reproducible shape
change and the rotated control explains nothing. So IgG's true motion IS partially
warp-like; flex_pca's state maps lose that signal. The headline finding is now
about the EMBEDDING, not about warps. `fromp/top` subsetting works fine
(`simple_parameters_phases.f90:655`); the raw difference is just dominated by a
Vbar-aligned term, which is real (a compact->extended hinge moves density out of
the core), not an artifact.

**Design fix from the literature.** In image-based LDDMM the momentum is
`m = lambda(x)*grad I(x)` — a SCALAR field times the gradient, zero in homogeneous
regions. So phi should be a scalar dictionary multiplying grad Vbar, not a 3-vector
one: cuts K by 3 and removes the degenerate directions by construction. And the
dictionary should be low-frequency elastic-network normal modes, not a B-spline
grid. See doc section 11.

**LINEARISATION VALIDATED (2026-08-05).** All these tests project onto the LINEAR
TANGENT SPAN at Vbar; DynaMight solves a nonlinear problem, so the obvious
challenge is that linearisation was the error. It was not. SimpleITK diffeomorphic
demons (fit on EVEN, scored on ODD) lands within **+0.003 / +0.014** of the tangent
projection (b00: 0.684 vs 0.681; b80: 0.565 vs 0.551) at **10-11 A mean
displacement — ABOVE lp/2 = 7.5 A** where linearisation should have failed. The
tangent methodology is sound. Script: `scripts/flex_pca_nonlinear_check.py`.
(A hand-written demons attempt DIVERGED first — scored below the identity at 27 A
displacement, no step normalisation or convergence check. Easy failure to repeat;
always check a registration beats the identity. A SITK B-spline FFD run returned
zero displacement = failed optimisation, not a result.)

Second finding there: both methods reach only +0.684/+0.565 against ceilings of
+0.979/+0.977, so **a properly fitted nonlinear diffeomorphism explains only ~half
the consensus->state difference** — independent support for the growth-term result
below, and proof it is not an artefact of linearising warp but not growth.

**WHY IT STILL WORKS IN RELION — the untested variants:** DynaMight moves ~10^4
GAUSSIAN PSEUDO-ATOMS, not a density field (a Gaussian centre is observable in all
3 directions; a density-field displacement only along grad Vbar — structurally a
different model); it fits to IMAGES with a jointly inferred latent; and it runs at
higher resolution where grad Vbar has internal structure, not just envelope.
**Scope of the negative: DENSITY-FIELD deformation as an encoder basis, at
lp=15 A, on IgG. Pseudo-atom deformation is untested and is what RELION's success
actually points at.**

**REOPENED as UNBALANCED TRANSPORT (2026-08-05, late).** Everything below tested
MASS-CONSERVING warps, spanning `{u.grad V}`. WFR / unbalanced OT (see the existing
`scripts/flex_models/unbalanced_transport.py` and the "unbalanced transport"
section of `python_first_flex_heterogeneity_models.md`) is explicitly
NON-mass-conserving and spans `{u.grad V + g*V}`. Adding growth columns `b_m*V`
alongside warp columns `-b_m d_d V`:

**The GROWTH term beats the WARP term at a THIRD of the rank** on both GT pairs
(growth(sp5) K=30: +0.245/+0.270 vs warp(sp5) K=90: +0.199/+0.222). rigid+growth
lifts rigid(D=3) from +0.149/+0.211 to +0.229/+0.261. On h00a, WFR(sp4) (+0.288)
and growth(sp4) (+0.281) **BEAT the linear bases** (+0.254 eigA, +0.248 bagged) —
the first generated basis to do so in absolute terms.

Physical reading: at lp=15 A in AVERAGED reconstructions a hinge motion is not a
clean displacement — the moving domain smears across positions, reading as density
fading out in one place and in in another, i.e. an occupancy/amplitude change.
Mass-conserving models cannot express that; WFR's growth term exists for it.

**What still blocks it is PARAMETER EFFICIENCY, not capability:** WFR needs K=252
where bagged linear needs K=5, and per-particle SNR is 0.0028, so K=252 is
unaffordable by ~2 orders of magnitude. Next thing to try: a LOW-RANK growth
parameterisation (few occupancy modes, not a dense bump grid). Untested.

**PROXY CAVEAT applying to every basis comparison in this whole investigation:**
they all score bases on explaining RECONSTRUCTED MAP DIFFERENCES, which is a proxy.
**The deformation encoder was never run end-to-end** (embed particles ->
weighted-backproject -> score states). `probe_external_basis` makes it cheap. Do
that before concluding anything final.

**Mass-conserving encoder case closed 2026-08-05.** Four generated parameterisations tried:
B-spline vector fields, scalar LDDMM momentum, ENM normal modes, rigid domains.
Against a PROPERLY HELD-OUT linear eigenbasis (basisrep basis A, fitted on eo=0
only, scored against the eo=1 reconstruction) the linear model wins at every rank:
eigA(5) gives +0.256/+0.254/+0.166 on the three GT bin pairs with **5** parameters,
vs rigid(D=3) +0.149/+0.211/+0.046 with **18** — ~3.6x more parameter-efficient AND
absolutely better.

Useful sub-result: **rigid domains beat ENM** (+0.211 vs +0.144 at matched rank on
the strongest contrast, a 46 % gain), and D=3 is the optimum matching IgG's
2-Fab+Fc architecture, with a rotated control near zero. So "generated basis" was
never the problem — a uniform elastic network on DENSITY-placed beads with a 39 A
cutoff was: it bridges the hinge gap and has no chain connectivity. **Anyone
revisiting this should start from rigid domains, not ENM.** Script:
`scripts/flex_pca_rigid_domain_basis.py`.

**The one thing still untested** is bias-variance: all the above measures BIAS on a
single linear fit, but nominally equivalent linear fits vary **5x** in explanatory
power (see [[flex-pca-basis-rank-over-specified]]) while a generated basis has zero
fit-to-fit variance. Scoring embeddings from N independent linear fits against the
fixed generated basis end-to-end is the only comparison that weighs both together —
and the cheaper alternative it suggests is simply AVERAGING several linear fits.

**Scope — narrowed on 2026-08-05 after the user pushed back.** What is falsified
is the LINEARIZED DENSITY-FIELD basis span{phi.grad Vbar}, as an explanation of
flex_pca's state-map differences. It does NOT falsify deformation modelling as
practiced by DynaMight/3DFlex, and must not be cited that way. Those fit the
deformation to PARTICLE IMAGES, and DynaMight in particular moves ~10^4 Gaussian
pseudo-atoms rather than warping a density field
(`relion/src/pipeline_jobs.cpp:5549`, "1-2 Gaussians per amino acid"). A Gaussian
centre is observable in all 3 directions; a density-field displacement is
observable only along grad Vbar. That aperture gap is the likely explanation for
the null, which makes a pseudo-atom parameterization a genuinely different and
untested proposal.

**Calibration (all with the same estimator).** synthetic hinge (known warp)
+0.76; Tomotwin-100 GT species pairs (known compositional, no embedding) -0.10 to
+0.09; Tomotwin-100 flex_pca states +0.07; IgG flex_pca states 0.00 to -0.06.
IgG sits with the compositional negatives. Tomotwin is a NEGATIVE calibration
point — it is 100 distinct PDB complexes, so a null there confirms the test
rather than invalidating it (the original doc reasoned this backwards).

**Two traps worth remembering.**
1. The eigenvolume subspace-angle test as first written is **underpowered and
   must not be used** — the masked band-limited space has only ~500 effective
   dof, so a rank-327 subspace captures ~80 % of anything and the phase-randomised
   null tracks the signal exactly. Score STATE DIFFERENCES cross-halfset instead;
   that makes noise-fitting score zero by construction.
2. The result scores flex_pca's OUTPUT, not IgG's true motion. Both halfsets
   share one embedding, one basis and one set of kernel weights — only the images
   differ — so a systematic embedding error reproduces at +0.947 exactly like real
   signal. That the states 4-5 difference is highly reproducible, is not a
   contrast change (corr with Vbar = 0.04), and resists any smooth rank-408 basis
   is itself suspicious and points back at [[flex-covariance-basis-defect]].

Instruments: `scripts/warp_state_difference.py` (the good test),
`warp_positive_control.py`, `warp_amplitude_sweep.py`,
`warp_subspace_angle.py` (the underpowered one, kept as a record). All run in
`.venv-flex-graph` in seconds against existing MRC outputs — no SIMPLE build.

Related: [[flex-pca-rename-and-cleanup]], [[igg-test-project-setup]].
