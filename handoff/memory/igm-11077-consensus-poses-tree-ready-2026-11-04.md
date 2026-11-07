---
name: igm-11077-consensus-poses-tree-ready-2026-11-04
description: "IgM tree input = the 683k consensus AS IS (user decision 2026-11-06 ~12:20; t2more_it3.simple, 3.85 A; baseline run tree/igm683k_s1) WITH the register caveat: ~1/3 mis-registered by 60 deg, J slot ~60 % filled, consensus locked on a heavy lobe ~150 deg from the J; register not recoverable for side views (RELION pilots v1/v2 failed); s6 = >40 A lopsidedness; regemd TRAP"
metadata:
  node_type: memory
  type: project
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-06T13:33:29.026Z
---

**DECISION 2026-11-06 ~12:20 (user): use the current consensus AS IS as the IgM tree input, with the register caveat.**
- Input: `tpick_all/3D_tests/T2/neigh683k/2_refine3D/t2more_it3.simple`, 683,255 particles, FSC 3.96 A (RELION 3.85 A).
  The earlier tree run `tree/igm683k_s1` (igm.simple = a copy, rows01.txt) used exactly this input and stands as the baseline.
- Caveat to state with every IgM tree result:
  - about 1/3 of the particles sit in the wrong 60-deg register (face-on 2D: 32 % off, mostly +-60);
  - the J slot is ~60 % filled (J/Fc .61 vs EMD .12-.24);
  - the consensus is locked onto a heavy lobe ~150 deg from the J, not onto the J chain;
  - a split that separates particles by J register, or a "J slot" state, is this nuisance, not biology;
  - side/tilted registers cannot be fixed by any fixed-pose method tried (RELION pilots v1/v2, every scorer);
  - hexamers are < ~1 %.
- Register work is closed. The paper-style multi-class refinement with re-alignment remains an optional longer run.

**PILOT v2 2026-11-06 12:10 (user-approved): + relion_particle_subtract (J-slot wedge, box 40) + K=3 seeded J/Fc/lobe from face-on
subset 1 -> FAILED for side/tilted views** (job 54988649): no class kept an empty slot (J seed drifted to Fc-like, CC .56 vs .71);
held-out face-on agreement .785 (= gap-only .783; all-j0 .678); tilted/side ties (median margin 0, >1 J copy .91-.98); called
registers 7.11 -> 8.84 A, 60-deg self-corr .34 -> .74. VERDICT: fixed-pose register reading works for face-on views only; the
paper's route (multi-class refinement WITH re-alignment) is the untried kind.

**PILOT 2026-11-06 11:00 (user-approved): RELION C6 symmetry expansion + one-slot K=2 classification FAILED** (`symexp/`, jobs
54987757 / 54987915; tb100k + 9,235 face-on with reference-free registers, poses in the dock3 frame, box 100).
- No signal subtraction. The classes split by the consensus HEAVY LOBE (~147 deg from the J, tied to the consensus frame, not to
  each J), not J vs Fc. Face-on agreement .008 (lobe-angle rule .386) vs all-j0 .678; called-register halves 7.11 -> 7.88 A and
  J/Fc .61 -> .64 (same signature as every earlier scorer).
- Validated on the way: SIMPLE->RELION export (CC .996/.973 vs rec_orig), the frame math, the 60-deg copy sense (.993). The
  face-on registered map has a clean J (J/Fc .20, five equal Fc). TRAPS: RELION 4 relion_preprocess cannot operate on STAR
  files (Fourier-crop in Python); RELION reorders output rows (map copies by image name + pose); the consensus "lightest slot"
  is not a usable J (three slots .65-.69), so take the J from the face-on registered map.

**VERDICT 2026-11-06 08:40: DJ r300k state 6 is NOT the registered pentamer; no de novo map so far carries the register.**
- Band check (job 54986056, `regfix/regband.py`, `reps/regemd/igm_regband.54986056.out`; core r < 130 A, CC per band vs EMD
  at k*60 and vs EMD-C6): at 60-30 A s6 k=0 .26 < C6-level s1 .30 (consensus .51, EMD-C6 .77), relative 60-deg drop .62 >
  EMD's own .37 (a noisy or register-mixed EMD copy cannot exceed .37) = lopsidedness that docking lays on EMD's J side. At
  30-20 A all 18 de novo class maps (ours + DJ's) are ~0 vs consensus .40: solve3D's "16 A" is shared-reference imprint.
- TRAP: regemd's 20 A-Gaussian contrast is a > 40 A readout. A register claim needs (i) a relative 60-deg drop <= EMD's,
  (ii) EMD-like CC at 30-20 A, (iii) docking into the consensus putting the J at the consensus J window (s6: untested; the 158-deg figure came from
  DJ's v2 docking, which failed its own gates, so it is not a measurement; the verdict rests on (i) and (ii)).
- `reps/djrep` replicates (54984558/9) lost their premise; cancelled by job ID on the user's instruction (2026-11-06 12:35),
  with rB's workers 54989115-24 (my own earlier scancel had been blocked by the auto-mode classifier). Next route (user's call): RELION C6 symmetry expansion + focused classification of one slot on
  the consensus poses, validated on the 9,235 reference-free face-on registers.

**CURRENT STATE 2026-11-06 03:25 (supersedes the class-based "breakthrough" of ~02:10):** class route NOT validated.
- EXCEPTION (05:45, DJ session): DJ random 300k (`dj_ws2/pool200/r300k`, box 200, 6-state solve3D, full mask, no 2D
  selection) state 6, 72,658 ptcls, 16 A: register contrast **+.138** (mirror -.105), a SHARP k=0 peak (core CC .740 vs
  <= .602 at k != 0), sectors J-first .31 1.03 1.29 .89 .59 1.20 (second dip at k=4: register mix?). First map above the C6
  band. s1-s5 (102k/34k/28k/34k/29k) -.008..+.028 (s3/s5 junk, core CC .36/.18).
  Half-map check (06:00, job 54983521, reconstruct3D nstates=6 state=6 in `r300k/2_s6_halves`): even +.138, odd +.137, k=0
  best in both (.739/.741 vs <= .604) = criterion met. BUT solve3D aligned both halves to ONE shared merged NU reference in
  every stage ([[refine3d-lp-set-scores-merged-map-trap]]), so this rules out reconstruction noise only, not reference imprint.
  07:35 independent check launched: two disjoint random 300k of the same DJ pool (rows 300,001-600k / 600,001-900k of the
  same seeded shuf), identical protocol + regemd both hands: `reps/djrep/rB` (54984558, csbdevel), `rC` (54984559, norm),
  ~9 h. A J-gapped state in either = s6 confirmed and the route scales (each draw = production chunk). Our own pool: 0/2
  draws (rep_A/B) gave one; DJ's: 1/1. DJ session tests s6 as the matched-reference scorer on tb100k (face-on agreement
  >= .85 vs gap-only .783, FSC not worse than 7.6 A) before any 683k fix.
  Band readout (08:28, pre-scoring checks of cancelled job 54986097, `ph_tests/s6score3`): placed on dock3 in EMD's register
  (contrast +.134 reproduced), s6's shell CC vs EMD .78/.43/.21/.09/.06/-.03 (120-60/60-40/40-30/30-25/25-20/20-16 A), vs
  rec_orig .58/.20/.07/.04/.05/-.01, yet its halves agree .95-.99 to 25 A and .82-.90 down to 9 A = shared-reference imprint;
  real resolution ~40-60 A, so +.138 is a >= 40 A feature. Scorer route held for the coordinator's band rule (job 54986056).
  Docking trap on the way: [[dock-volpair-wrong-optimum-check-copies]].
- rep_A/rep_B (disjoint random 300k of the 2.4M, solve3D 6 states, full mask): A5 145k / B2 117k looked J-gapped (.25-.28)
  only because regfix/jslot_self.py took the disc normal from second moments of the WHOLE thresholded map (arms dominate;
  normals 45-50 deg apart for maps docked into one frame). Core-only normal (r < 125 A, regfix/jdir2.py): A5 .38, B2 .19,
  both as GRADIENTS (off-centre-ring pattern), not EMD's single empty sector (.12 1.03 .93 .89 .95 1.20).
- Deposited-frame test (reps/regemd, regfix/regemd.py): dock into EMD-13921 through a 280 A core mask at 20 A, read on EMD's
  axis/centre; register contrast = core CC(k=0) - best CC(k!=0). EMD .189, EMD J-filled .075, EMD-C6 .001; consensus -.020,
  A5 +.013 (mirror +.005), B2 -.061 (mirror +.032), all 12 replicate classes -.034..+.038 = C6 level. Caveat: our cores match
  EMD worse than EMD-C6 does (CC .73-.79 vs .895), so "register-averaged" vs "different near the core" is not separated.
  Full-mask docking into A5: EMD-C6 fits better than EMD (.864 vs .818).
- DJ synthetic: with a MATCHED reference the register is recoverable per particle (.95; side views as informative as
  face-on). The obstacle is the reference, not the information. EMD-13921 is from this dataset's paper (Chen 2022: 5 ab
  initio classes + heterogeneous refinement, 2 kept).
- 2D candidates (fullreg/states.txt, 1,643,183) keep 90 / 94 % of A5 / B2 and 8-21 % of junk-like classes (reps/clsoverlap).
- TRAP: any J-sector readout needs a core-only normal AND a C6 control read in the same frame; centring alone makes gradients.
  The |F| n10 test inherits the flaw: tree/symfold.py also takes the normal from the whole thresholded map.

**UPDATE 2026-11-06 ~00:20 (DJ session; adds to the state below):**
- Gap-only reference G (consensus C6 Fc ring minus the J wedge, `ph_tests/gaponly/G_gaponly.mrc`).
  - Matches the reference-free face-on registers .78 (.91 at margin >= median).
  - On tb100k it only partly empties the J slot (.59 -> .46) and pulls the map toward C6. Use it as a face-on/confident scorer, never as a refinement reference.
- Six-lobed face-on 2D classes (11/14/18/19, 7.7k) are NOT hexamer evidence (`ph_tests/sixlobe/`).
  - >= ~75 % carry a J gap as visible as the gapped classes'. Cross-fitted J/Fc .40 vs .21 (all), .27 vs .11 (confident); no gap = 1.
  - <= ~25 % show no gap: hexamers or gap-invisible pentamers, not separable by this test.
- DJ class-based de novo (51 side-view classes, 410k, `dj_ws2/checks_3D_dj/`): 3 states (187k/88k/135k, 16 A), all at the C6 level.
  - Deposited-frame register contrast (regemd, 03:28): s1 +.033, s2 +.034, s3 -.042 (mirrors +.022/+.012/-.003); core CC vs EMD .63-.77.
  - solve3D_addon (+292k from 48 more classes -> 264k/211k/227k, 04:00): still C6 level, +.034/+.029/-.044 (mirrors +.025/+.009/-.007).
  - The 01:35 "lean pentamer" call (jslot .36/.42, |F| n10 .26-.45) is WITHDRAWN: both readouts took the normal from the whole map.
- TRAP: 2D slot windows taken from an average's own six-fold phase (face_register.slot_profile) have a sign error and drift between averages; the face-on profile has 12-fold sub-lobes.
  - Use fixed windows at the 3D slots (17.75 + 60k, J 197.75 in the dock3 frame). The stored 3D windows are verified: EMD as-is .13, consensus .60.
- TRAP: consensus-pose averages can carry register noise fitting. Cross-fit instead: call the register from one Fourier half (kx+ky parity) and average the other. The best-of-six placement does not depend on the consensus register.

**CURRENT STATE 2026-11-06 ~00:00 (supersedes the verdicts below):** user: the 3.85 A map "looks like a hexamer". Slot masses
about the axis: J slot / Fc slot = .59-.70 in our consensus vs ~.11 in EMD-13921 (9 vs 70-91) -> the J chain is smeared over
the ring = core register scrambled (or hexamers). The tb3 "97.5 % keep register" and real7 "no hexamers" tests scored
particles against references derived from the consensus itself: self-referential, they cannot detect a register average
baked into the reference (the lopsided rod/Fab mass sets the vote). The register-fix blur (7.4 -> 12.3 A) does not prove
registration (a 1 px axis error blurs too). Solid: 2D face-on (reference-free) = 9.4k particles with a J-chain gap (pentamers);
averaged with their CONSENSUS poses they lose the gap (n10 .25-.40 -> .04) -> the consensus scrambles their register.
Decisive route: reference WITH a clean J gap, correct hand (EMD as-is fits better), J/gap/rod mask; success = J slot empties
in independent half-map reconstructions. Never call the consensus "registered" from self-referential scoring.

Register methods tried 2026-11-06 (all on tb100k): full-mask refine=neigh from the consensus (automsk=no, lp 10) STAYS at J/Fc
.60 and drops 7.4 -> 10.2 A; refinement from EMD with the J slot filled: 24 A, J/Fc .87; per-particle scorers with EMD (either
hand; its Fab-1 rod-root heavy slot sits at J-60, ours at J+120) or a gap-only reference scramble side/tilted views. Only
face-on gaps are readable (reference-free, ~4 % of IgM). Next route = class-based selection (multi-state full-mask de
novo, keep classes with a clear J gap; regfix/jslot_self.py: EMD .13, C6 1.00, consensus .63).

**RETRACTION 2026-11-05 ~23:00 (itself superseded above):** the "6-fold register average" verdict below is WRONG. About the docked pseudo-6-fold axis the
tb100k consensus map correlates with itself rotated 60 deg at only .34 (30-8 A) / .04 (12-7 A); after a wrong-hand register
"fix" scrambled it, .70 / .21. So the consensus is mostly in one register; its 6/12-fold |F| and ring n6 come from the
pseudo-hexagonal core (5 Fc + J chain). The per-particle register fix with a docked EMD-13921 (mirror hand) made it WORSE
(7.4 -> 12.3 A); EMD-13921 as-is fits the consensus better than the mirror. A uniform rotation of all poses can NOT test the
axis (it just rotates the map). Input-free retry: score against the consensus's own rotations, cross-validated halves.
RESULT (regfix/tb3, 23:00): scoring tb100k against the six 60-deg rotations of the consensus itself, eo=0 vs the odd half and
eo=1 vs the even half: j=0..5 = 97549/937/317/160/204/833; margin >= median: 49964 j=0 vs 38 other; FSC 7.43 A orig, 7.60 A
all calls, 7.43 A confident only; self-corr unchanged. => the consensus is registered; no register fix needed. Pentamer vs
hexamer stays open (de novo multi-state cannot separate them; per-particle 7-way scoring with a matched hexamer is the test).
Axis checks (23:30): moving the axis point 2 px perpendicular (4 directions) RAISES keep-register to 99.4-99.5 % (0-3
confident switches) - the dock3 point gives the most switches, so it is the best point tested; axis-independent: the max
60-deg self-corr over an axis search is .44 (rec_orig) vs .82 scrambled. Verdict robust: >= ~97.5 % registered.
PENTAMER vs HEXAMER (23:45, DJ session, ph_tests/real7/): cross-validated 7-way scoring of tb100k (six rotations of the other
half map + a hexamer built from that half) -> D (H minus best register) median -1.49, no H-like mode; called H 1.2 % (D > 0),
.22 % (D > .33) = at/below the synthetic pure-pentamer false-H rates (.028/.002). The 2,000 most H-like particles
reconstruct to the same lopsided pentamer as the 2,000 strongest P (J/Fc .70 vs .69, 60-deg self-corr .565 vs .569, halves
agree). => hexamers < ~1 % of the consensus set. Not covered: particles outside the consensus (discs, dropped classes).
Original (wrong) text follows.


`/mnt/beegfs/elmlund/afan/igm_cryolo/tpick_all/3D_tests/T2/neigh683k/2_refine3D/t2more_it3.simple` holds 683,255 particles.
- The poses come from refine3D refine=neigh with automsk=yes, which is a core-only mask (radius 140 A or less), on top of a de novo solve3D from the user's 15 core 2D classes.
- 12 of those classes are side views.
- FSC 0.143 is 3.96 A, and the RELION reconstruction of the high-res iteration-4 poses gives 3.85 A.
- The map is **6-fold**, measured as ring-centred angular power down the disc normal at 110 A:
  - our map: n6 0.61-0.79, n5 0.08-0.23;
  - EMD-13921 (deposited pentamer): n5 0.76, n6 0.04, and it shows the asymmetric subunit-1 F(ab')2 rod held by the J chain.
- The pentamer's Fc core is pseudo-hexagonal, so a single model averages pentamers over 60-deg registers.
- The old blind control (greedy re-search) only tested 72/144-deg flips. Greedy also stays near the existing pose, so it cannot detect this error.

Test bed: 100k random particles at box 200 (`tree/tb100k`), de novo, full mask, automsk=no.
- 1 state: drifts back to 6-fold (n6 0.73).
- 3 states (`tree/tb3d_s3`):
  - 43.7k particles: no ring symmetry (junk or side views);
  - 29.9k: 5-fold with the rod by ring centring (unconfirmed by |F|);
  - 26.4k: 5-fold by ring centring (unconfirmed by |F|);
  - no 6-fold class.
- The paper (Chen et al. 2022 Nat Commun, PMC9588798) says hexamers are under 5 % and that its SEC step removed them. Their recipe was 5 ab initio classes plus heterogeneous refinement, keeping the 2 best classes.

- 2026-11-05 face-on 2D (DJ session, `ph_tests/`): fine 2D of our 27,039 face-on stars gives 9.4k pentamer classes (J-chain gap,
  |F| n10 .25-.33, Fab-1 rod at the gap) and 8k six-lobed classes; C6-symmetrised, six-lobed vs pentamer .983 vs a synthetic
  hexamer's .90. CORRECTED same night: that synthetic hexamer was built about a 5.4 px-off axis; built on the true axis it is
  .965 vs C6(pentamer) (noise-free), so the real .983 only LEANS to register-averaged pentamers (suggestive, not decisive).
  2D CAN resolve the register (9.4k gapped pentamer face-on particles).

**Why:** the "tree-ready" verdict was wrong. Every tree result on these poses (root junk split, r1 leaf, core-only sites) reflects the register average.
**How to apply:**
- For IgM, separate classes with a full mask and several states before refining.
- Never trust a single-model IgM consensus without the register caveat (the user accepted the 683k consensus as the tree input on
  2026-11-06 with that caveat; see the DECISION section at the top).
- Check 5- vs 6-fold with the CENTRING-FREE test (`audit/fold_ft.py`: angular power of the 2D Fourier amplitude of the
  disc-normal projection; EMD-13921 = n10 .32 vs n6 .01 at 45 A). Ring n5/n6 at 110 A flips with the chosen centre (EMD-13921:
  .05/.14 at COM, .74/.01 at max n5+n6), so it is unreliable alone; by |F| the consensus stays 6-fold, but the tb3d_s3
  "5-fold" classes show no clean n10 (unconfirmed, 2026-11-05).
- The 2.4M box-200 pool is `pool200/igm_pool_sc.simple`.
- Related: [[igm-retain-side-views]], [[box600-simple-master-needs-64g]].
