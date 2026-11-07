---
name: impsamp-toy-campaign
description: "2026-10-04 flex2D- vs abinitio2D-initialised 3D importance sampling on toy is1; flex2D wins via junk removal; greedy 2D-score ranking beats random-within; paths, scripts, results"
metadata:
  node_type: memory
  type: project
  originSessionId: 95245f2d-90ec-4676-978e-28a5533c3f54
  modified: 2026-10-04T11:54:36.638Z
---

User ask (2026-10-04): test flex2D-initialised vs abinitio2D-initialised importance sampling for 3D, optionally with weights
derived from the flex2D embedding; develop on a toy first, then demonstrate on a molecule.

**Mechanism (code):** abinitio3D (balance=yes default) draws each iteration's subset with sample4update_class
(src/main/ori/simple_oris_sampling.f90): equal share per ptcl2D class over cls2D-selected classes, members ranked by ptcl2D
'corr' (greedy top-n in early stages; frac_best 0.5 / 0.85 later). Only ptcl2D class / corr / state / shifts reach 3D;
class averages never do. NSAMPLE default 10k per state (FULL_SAMPLE_SWITCH 0.9: small sets get no subsampling at all).

**Toy is1** (batch2 /mnt/beegfs/elmlund/afan/impsamp/is1): 60k ptcls, 50S EMD-8440/8445 80/20, half within a 15 deg cone,
10 % pure-noise junk, SNR 0.05, box 160 @ 2.62. Particles are stored ORDERED (state1, state2, junk) -> a worker range can be all
deselected -> abinitio3D worker THROW "no active particles to sample" and the master hangs (flagged as task; use NPARTS=9 there).
2D: a2d60 (seed), a2d180 (control), flex180 (seeded, deselects 5,980/6,000 junk). Scripts in scripts/flex2D/impsamp/ and on
batch2 impsamp/: is_gen.sh, is_a2d.sh, is_flex.sh, is_ab3d.sh (NPARTS env), is_relabel.sh (ptcl2D class/corr/state swap via
import donor + replace_project_field + selection infile, verified 0 mismatches), is_sampler_sim.py, is_groups_proto.py, is_make_labels.py.

**Results (2-state abinitio3D, purity / rare recall / junk active):** flex default B 0.965,0.956 / 0.75,0.72 / 0.3 %;
a2d default A 0.821,0.825 / 0.17,0.19 / 100 %; flex random-within C 0.927/0.57; embedding HDBSCAN groups random-within E
0.927/0.57; uniform on flex project G 0.907,0.970 / 0.46,0.77 (high variance); a2d random-within D 0.800/0.00 (junk took a state).
Python composition prototype said greedy corr ranking cuts the rare share (a2d 0.07, flex 0.13 vs 0.20), yet greedy B beat
random-within C/E in 3D: clean early subsets matter more than rare share. Latent kNN is strongly state-homophilous (0.84 vs 0.007)
but inverse-density / HDBSCAN weights did not up-weight the rare state (oracle 0.44). Round 2 running: H = a2d with junk classes
deselected (isolates junk vs partition), E' = embedding groups with flex corr ranking; F oracle pending.
Related: [[judge-2d-partitions-by-heterogeneity]], [[prototype-in-python-before-fortran]], [[flex2d-findings-2026-09-21]].

**Real data (2026-10-04 evening).** EMPIAR-10076 fair pair, 3 runs each (4-state abinitio3D, snapshot 6c8ca8f15, identical
65,949 flex2D-pruned particles, published-label purity): flex2D 1000-class partition 0.710 / 0.692 / 0.704 vs abinitio2D 1000
0.503 / 0.486 / 0.577 -> no overlap. Caveat: the particle set was pruned by flex2D evidence (bottom 50 %, 75 % of unlabelled
removed); pure abinitio2D 200 on all particles gave 0.418. TmrAB (EMPIAR-10317, 1,547,017 ptcls, labels 1/2/4 view-independent,
white protein; impsamp/tmrab, 200k random active): a2d100 + a2d300 done (4 h each, slow 405 GB stack), flex300 running.
**Purity ideas tested and negative:** latent k-NN smoothing of 3D labels lowers purity (toy 0.965 -> <=0.96, 10076 0.710 ->
<=0.709); flex2D class vote 0.68. score_states (flexref5 build, fixed pose + shift search) on ab3d_prune1000_r2: posterior
argmax == current label 99.96 % self-scored and 100.00 % cross-half (even maps -> odd ptcls and vice versa) -> purity unchanged
0.692. Fixed-pose scoring cannot reassign states when each particle's pose was refined against its own state map; the four
abinitio3D state maps correlate only 0.52-0.72 unaligned at 20 A (not clearly co-registered). A valid test needs per-state
orientation search (pose marginalisation), not just shifts. Scripts: scripts/flex2D/impsamp/e76_score_states.sh, e76_score_xhalf.sh.

**TmrAB flex2D blur diagnosis (2026-10-05 night, user asleep: "see if you can find the issue").** flex300 averages are visibly
blurrier than a2d300 (a2d shows TM helices at 6-9 A; flex 8-16 A); 2D state purity is equally poor (flex 0.479 / NMI 0.033, a2d
0.494 / 0.039, chance 0.36). Coarse class view/in-plane spread vs deposited poses is the same for both (~43 deg). flex300's polish
phase is INERT: rotation changes >2 steps 0.000-0.001, mean shift increments 0.02-0.35 px in every polish iteration (stages 7-9),
i.e. poses are essentially the 20 A label-phase poses. Suspected mechanism: gold-standard half templates (~340 ptcls/half) +
apply_template_filter (2f/(1+f) of the HONEST template eo-FRC, simple_flex2D_model.f90) strip the high-res template content, so
polish has no high-res alignment signal; cluster2D's merged references inflate f and keep it. Variants launched from a2d100
(impsamp/tmrab, is_flex_var.sh): flex300_nofreeze (class search in polish), flex300_gold0 (FLEX2D_GOLD=0), flex300_fm2
(FINAL_MERGED=2), flex300_notf (FLEX2D_TMPL_FILTER=0). Referee map: refmap/tmrab_ref100k.mrc from deposited star poses of 100k
particles outside both 2D subsets (shift sign chosen by 20k FSC).
Map-based FRC vs refmap/tmrab_ref100k.mrc (deposited poses, 100k ptcls outside both subsets, FSC 8.9/8.6 A): a2d300 30-12 A 0.520,
12-8 0.079, res@0.5 19.7 A, 21 % of ptcls in classes reaching 0.5 by 12 A; flex300 0.375 / 0.056 / 23.0 A / 0.7 % -> the deficit is
real and already at 30-12 A (alignment precision, not just fine detail). cluster2D_distr refine=inpl (simple_private_exec; refs=
flex cavgs, flex sigma2_state.bin copied; cls_init=prev is rejected) was INVALID: refs= does not keep the project's labels; classes
came out scrambled (purity 0.363 = chance, NMI 0.001, map FRC 30-12 0.20) while poses did not move (in-plane dist 0, shift 0,
converged after 1 iteration). The 'self-consistent fixed point' reading rests only on flex2D's inert polish statistics. abinitio2D escapes by co-evolving classes, poses and references while
marching 15 -> 8 A with stochastic search. Cancelled F3/A3 and the fm2 variant (QOS cap; diagnosis first).
Trap: a run.sh that `touch DONE` even on failure leaves a stale DONE; wait on RC= in the log, not on DONE.
2026-10-05 morning: the four TmrAB 3D arms (ab3d_F1/F2/A1/A2) TIMED OUT at is_ab3d.sh's 8 h limit (iteration ~42, no score) under the
QOS-cap contention, and their already-submitted simple_tmrab workers kept running as orphans; user: "kill the jobs that aren't useful,
using a lot of cores for no reason". Cancelled them plus the four TmrAB 2D diagnostics (nofreeze/gold0/notf/a2dpolish) in favour of a
fast toy. Lesson: give multi-hour masters a generous -t, and when a master dies, cancel its workers. flex2D's polish is inert on
RNase P (labels 22 A) and the spliceosome (44 A) too, not only TmrAB: poses are always the label-phase poses; small or featureless
particles (TmrAB in its belt, hemoglobin) pay for it. Toy for the fix: impsamp/pt1 (30k ptcls from tmrab_ref100k at box 128 /
2.154 A, uniform views, SNR 0.02, shift +-5 px, mask 180 -> 20 A label band); evaluation pt_eval.sh (exact-map FRC + in-plane
consistency within 5-deg true-view neighbourhoods + polish search stats).
**TmrAB fix campaign on tmf (2026-10-05 afternoon).** tmf = the 30k tmini particles copied into one contiguous stack
(impsamp/tmf/stack.mrcs; seed 20 min vs 96 min, flex60 ~8 min); truth.txt = deposited poses/classes; eval = pt_eval.sh with
MAP=tmrab/refmap/tmrab_ref100k.mrc, EVAL=<subdir> (beware: a mid-line # comment once commented out A=/mkdir/cd). Map FRC 30-12 / 12-8 A:
a2d60 0.507/0.069; flex60 0.343/0.044; nofreeze 0.336; GOLD=0 0.382/0.057; FINAL_MERGED=2 0.374; TMPL_FILTER=0 0.344; SEED_EXPLORE
0.328 (hurts); HPSH=40 0.312 (hurts); ENVMASK 0.259/0.001 (polish moves 12 % of rotations but wrongly); LPLABEL 15 0.400, 12 0.417,
12+GOLD0 0.436/0.070, **10+GOLD0 0.456/0.072 (best flex)**; flex30 0.414 vs a2d30 0.542 (population helps both); hybrid abinitio2D
cls_init=prev from flex60 0.460/0.090 but reshuffles (19 % keep their class, NMI 0.21). Weighted averages (tmini_wavg.py,
W1/W2/drop) = no gain (weights ~1). Re-averaging gains: flex 0.343->0.400, a2d 0.507->0.531 -> ~0.02 generic averager loss + ~0.03
flex-specific (suspect flex's sigma2 in make_cavgs ml_reg). flex in-plane angles sit on a 1.67 deg grid (216 values) vs a2d
continuous. Polish search is global (sweep all rotations + shift refine) yet inert: the label-band optimum is a fixed point.
Deposited-pose in-plane spread does not discriminate (3-5 deg everywhere). All TmrAB partitions separate the deposited conformations
equally poorly (purity 0.44-0.46, chance 0.36).
**tmf follow-up (2026-10-06):** LPLABEL 8 + GOLD0 (60 cls) 0.470/0.080 (beats lpl10g0 0.456/0.072); LPLABEL 10 + GOLD0 at 30 cls
0.526/0.103 vs a2d30 0.542/0.064 -> at matched class count flex now wins 12-8 A and is within 0.016 at 30-12 A. Finer label band +
merged templates is THE TmrAB fix. **Spliceosome sp3d:** independent 3-state abinitio3D (flex_r1, a2d_r1) stalled at ~20 A FSC per
state (a2d_r1 cancelled); the stage lp schedule comes from the selected classes' frc2D (set_lplims_from_frcs) so flex got coarser
stages (19.4/15.9/9.8 vs 14.3/9.7/8.6) -> faster, unequal. Docked reruns (flex_dk, a2d_dk: multivol_mode=docked force_lp_range=yes
lpstart=20 lpstop=8, via AB3D_EXTRA in sp3d_run.sh) BOTH crashed at the split (stage 6): "insufficient previously sampled
class-balanced update candidates" (simple_oris_sampling.f90 sample4update_class sampled_only; still present on master); split cohort
29252/29276 trimmed vs per-class allocation. Abinitio3D in snapshot d1a20eff6 has no balance switch (class-balanced always on).
