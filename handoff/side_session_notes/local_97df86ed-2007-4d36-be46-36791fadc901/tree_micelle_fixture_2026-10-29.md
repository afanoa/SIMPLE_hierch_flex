# PfCRT-in-a-micelle fixture for hetero_tree (session tree-micelle-synth, 2026-10-29)

Everything is on beegfs: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/pfcrt_micelle/` (W below). Nothing is committed.

## Fixture (reusable)
- **Real inputs, read only.**
  - `pfcrt_core/3_refine3D/PfCRT.simple` (56,809 particles, core-aligned to EMD-20806) is copied to `W/real/real.simple`.
  - `W/real/ptcl3D.txt` and `stk.txt` are dumps of it. The stacks are dibyajyoti's 6_extract, box 240 at 1.3 A, ctf=yes, **Cs 0.001**.
  - `W/real/1_reconstruct3D` is an UNMASKED consensus from those poses. Its FSC is ~0.82 at 9.2 A and 0.03 at 7.8 A: the lp-8 pose limit.
- **Model** (`W/scripts/micelle_fixture.py`; `mrcio.py` holds the I/O):
  - **Protein:** 6UKJ atoms (PfCRT plus Fab Fv; the real particles show no Fab constant domains) as a Gaussian-atom map, B 60.
    - It is scaled to the real consensus by regression on the extramembrane protein (|h| > 28 A), with no offset: scale 6.45.
    - Membrane frame: normal = mean TM-helix axis, centre = TM centroid.
  - **Protein states:** the half of the transporter on the e1 < 0 side rocks about the in-membrane axis e2 through the TM centre.
    - The moving half is 6 TM helices, whole-helix assignment: 89-112, 126-148, 210-237, 241-259, 279-291, 314-340.
    - The Fv sits on the static side.
    - Rocking by 8 / 16 deg moves the half by a mean of 2.4 / 4.8 A (max 5.1 / 10.2 A).
  - **Micelle:** the REAL belt. It is the consensus minus the scaled protein, at <= 10 A, outside the union of the state masks, in the slab |h| < 38 A.
    - The belt is 1.5x the protein mass.
    - 256 bins each warp it continuously:
      - outer radius scale (inner radius 20 A fixed), sd 0.18;
      - ellipticity, sd 0.12;
      - one-sided thickening, sd 0.12;
      - thickness, sd 0.09;
      - tilt, 6 deg per in-plane axis;
      - density, sd 0.15.
    - These are MICAMP 1.5 x the SD table in the script. Micelle mass varies about ±35 % across bins.
    - The rest of the non-protein density (above and below the slab) is a static extra.
  - **Images:** simulate_particles of the protein state plus simulate_particles of the micelle bin, noise-free, then summed (the simulator is linear).
    - Real poses (shifts halved) and real CTF.
    - Then the OLE-null quality model: B ~ U(0,200), amplitude lognormal(qsig), white noise sigma.
    - Noise, amplitude and qsig were MEASURED on 2,000 real particles cropped to box 120 / 2.6 A. Results (`measure.json`):
      - a_sd_rel 0.18, qsig 0.29, padded SNR0 0.00117;
      - real-vs-model band corr 0.15 / 0.10 / 0.07 / 0.07 at 120-30 / 30-20 / 20-12 / 12-8 A.
  - **Datasets** (box 120, 2.6 A, 56,808 particles). They share every particle's micelle bin and quality:
    - `W/micnull`: one protein state.
    - `W/mic3`: 0 / 8 / 16 deg at 45 / 35 / 20 %.
    - Both: refine3D (prob, lp 8, 4 its, mskdiam 200) from 10 deg / 2 px perturbed truth -> `<ds>/src.simple`.
    - Truth: `<ds>/truth_full.txt` (row state micbin bfac snrq lam e2 e1 nu tilt dens micpc1-3); the full bin parameters are in `<ds>/mic_bins.txt`.
  - **Launcher:** `W/scripts/run_micelle_sim.sbatch` (env DS, STATES and POPS colon-separated, MICAMP, NB, SEED). One job per dataset, ~1-2 h.
- **Tree launcher:** `W/scripts/mic_gate.sbatch` is sst_gate4 that also records SIMPLE_TREE_* in split_env.txt and run_manifest.txt.
  - PfCRT keys: `nparts=5 pgrp=c1 mskdiam=200 minpop=3000 site_layer=yes site_states=yes`, ARMKEYS `kway_maxdepth=2`.
- **Noise.** White noise at the measured total variance had the wrong spectrum. Its band correlations were 0.22 / 0.14 / 0.05 / 0.03 against the real 0.15 / 0.10 / 0.07 / 0.07.
  - The tree inputs use band-calibrated colored noise instead (`micelle_fixture.py calib` / `merge_cn` -> `<ds>/all_cn.mrcs`).
  - Synthetic band correlations are then 0.149 / 0.100 / 0.068 / 0.064 against real 0.151 / 0.099 / 0.069 / 0.069.
- **Input poses.**
  - refine3D cannot hold poses at this SNR. From 10 deg / 2 px perturbed truth the median error grew: prob, white noise 41-47 deg; prob, colored noise 113-115 deg; neigh with the protein mask 106 deg (neigh scanned 100 % of the space).
  - There is no common rotation, and true poses reconstruct to Nyquist (corr 0.997 with truth), so this is not a convention bug.
  - Tree inputs are therefore the truth perturbed by oriops angerr 3 / sherr 0.75 (median 1.41 deg, 0.92 px), with no refinement: `<ds>/src_in.simple`, poses `W/poses/import_a3_s0.75.txt`.
  - That perturbation was chosen so the input consensus reaches FSC 0.5 / 0.143 at 8.9 / 7.6 A, matching the real PfCRT input (8.0 / 7.8 A).
- **Scorer:** `W/scripts/micelle_score.py RUN FIXTURE [--maps]`; `W/scripts/score_trees.sbatch` (env RUNS) adds the root pose error by step.
  - Per split node: children, state shares, one-vs-rest AUC per micelle parameter, the state_split re-test and guard lines.
  - Per leaf: composition.
  - NMI of the leaves vs state, micelle PC1 and quality, and the null bar.
  - With --maps:
    - each leaf map's fitted rocking angle and micelle scale;
    - the stage-B difference energy in the truth belt vs protein region.

## Diagnosis: 6 root-only trees (maxdepth 0, production split env, seed 1; W/trees/<run>/score.txt)
Arms: recipe17d no mask | recipe17d_micsup + truth protein mask as tree_mskfile | recipe17d_micsup + SIMPLE_TREE_SUPPORT=order.
- **The tree's own re-alignment scatters the poses.**
  - The root's median / p90 pose error goes 1.4 / 2.1 deg at input -> 8.2-8.7 / 137-140 after stA0 -> 9-18 / 143-153 after stage A -> 7.0-8.9 / 129-145 after stage B.
  - That is more than 10 % of particles flipped by local prob_neigh at lp 15/15/10, in every arm.
- **Every split is aligned vs misaligned particles.**
  - Stage-B groups differ in pose error: small groups at a median of 21-60 deg vs about 6 deg; poserr AUC 0.59-0.76.
  - Group maps differ by |diff|/|cons| 0.4-0.95, reproducible across halves, re-test z 15-73 vs z_crit 3.01.
  - The micelle never drives a split: every micelle-parameter AUC (radius, ellipticity, thickening, thickness, tilt, density, PCs) is 0.48-0.52.
  - Protein states are not separated either (mic3 NMI <= 0.004).
- **No mask:** the difference is diffuse (77 % of its energy outside both the truth belt and the protein; guard L/L_ref 0.26), so the guard merges it.
  - The null passes, but by accident. mic3 also stays one leaf.
- **With a mask, user or data-derived:** the difference sits inside the mask (82-99 % in the protein region), so the guard keeps it (L/L_ref 0.87-1.03).
  - The nulls fail: micnull umask K* 3, largest leaf 58 %; order K* 2, 54 %.
  - mic3 splits K* 2-3 with no state separation.
- **The order support itself works as a support:** on micnull it keeps 100 % of the truth protein mask.
  - Belt kept by radius 20/30/40/50 A: 0.95/0.49/0.01/0.00. Real PfCRT halves: 0.90/0.39/0.02/0.00.
  - It does not fix the tree, because the failure is the poses.
- **Real PfCRT parallel:** tree_umask2_s1 (another session's run: 6UKJ mask, rho floor 0.03) split its root hard (z 15.1, L/L_ref 1.08) and made 20 leaves; tree_fixpose_s1 stayed one leaf. That fits the same mechanism.

## Fixed poses: 4 roots (recipe17d_micsup_fp = micsup + recipe17b_fixpose/fixposes.diff; fixposes=yes)
- **Inputs** `<ds>/src_in_fp.simple`: as src_in, plus a balanced eo split (seed 41) and proj=1.
  - refine=cont refuses imported poses that lack them; only cont reads proj.
  - The split test's view strata come from the Euler angles.
- **Poses** stay at 1.41 / 2.06 deg (median / p90) at every stage.
- **micnull nomask:** re-test 0.047, z 2.52 < 3.01 -> K* 1, PASS; stage-B groups random.
- **micnull order:** cross-half EM collapses to one group -> K* 1, PASS.
  - The re-test line prints `****` and z Infinity: cosmetic, there is no pair.
- **mic3 nomask:** z 3.91, K* 2 (30066/26742), guard 0.999 kept. **A micelle-size split:**
  - radius AUC 0.58 / 0.42, micelle PC1 0.58 / 0.41, leaf radius scale 1.026 vs 0.976;
  - protein state V 0.011.
  - In micnull the same micelle sits just under the bar (z 2.52).
- **mic3 order:** z -0.65 -> K* 1. The order support removes the micelle split.
- **Protein states 0 / 8 / 16 deg:** separated in no arm (stage-B NMI <= 0.002).

## Oracle split tests (state_split with stage-B states = TRUTH labels; fixed input poses; W/oracle/<set>/split_<arm>)
- The first attempt was void: import_particles resets state to 1. Labels now go in via `prg=selection infile=labels`.
- **mic3 0 / 8 / 16 deg:** truth start board pairs 0.52 / 0.63 / 0.79 vs shuffle threshold 0.096, so the states are distinguishable at 39-10 A.
  - The cross-half EM (<= 19.5 A) changes 53 % of labels; the 8-deg group goes to deviation FSC 0.
  - Result: degenerate K* 2 (54760/2048 nomask, 55343/1465 order).
- **mic3 0 vs 16 deg:** EM changes 28 %.
  - Nomask focus: z 2.26, FAIL.
  - **Order focus: z 4.70, PASS, guard 1.21, K* 2 = the truth groups exactly (25643/11416).**
- **micnull micelle-radius halves:** EM changes 32 %.
  - Nomask: z 3.05 (marginal pass). Order: z 11.97.
  - Guard 0.575 in both, kept: the guard does not stop micelle splits.
  - The order support keeps the inner belt (95 % of 20-30 A, 49 % of 30-40 A), so a seeded micelle split passes more strongly under it.
- **Read:**
  - The tree misses the protein states because the cross-half EM band (<= 19.5 A) is dominated by micelle and quality and relabels away from the truth; a mask does not change that.
  - The order support helps a real protein difference pass the re-test, but is no protection against micelle splits.

## Pose-stability survey (read only; W/flips/flips.tsv, flip_survey.py)
Fraction of a node's refined rows whose pose turned > 30 deg from the node input, at stA0:
- qtoy2 <= 0.0007 (19 nodes); IgM <= 0.014 (78); OLE s1 / s2 <= 0.0245 / 0.0239 (roots; 33 / 34 nodes); qtoy null 0.0009; OLE null 0.0033.
- RNase P: no stA0 project survives, so the stage-B value is an upper bound:
  - 692k f692_r17d_s1: root 0.094, every other node <= 0.037;
  - r1211 <= 0.0007; r111A <= 0.0015.
- Fixture roots 0.289-0.314.
- **Real PfCRT tree_umask_s1 / tree_umask2_s1: root 0.407** (median turn 13 deg; stage B 0.53-0.56); depth-1 nodes 0.06-0.14.
- Gate threshold X = sqrt(0.094 x 0.289) = 0.165: the geometric centre of the root gap, a 1.75x margin on both sides.

## Remedies (both default off; recipe17d untouched; source in worktree tree-micelle-synth-24763a, uncommitted)
- **SIMPLE_TREE_SUPPORT=order** (order_support in simple_hetero_tree.f90):
  - The root mean's halves give cross-half power in the finest third of the FSC >= 0.5 shells (Gaussian 6 A window); Otsu on log power over the envelope; components >= 5 % of the largest; 2 A dilation plus 3 A cosine edge.
  - The result is node_r/support_order.mrc, used as tree_mskfile.
  - Builds: recipe17d_micsup (+ tree_mskfile pair), _fp (+ fixposes), _fpg (+ pose gate). All passed the fast gate (15/15).
- **SIMPLE_TREE_POSEGATE=X** (pose_flips + S_STA0 hook):
  - At the ROOT's stA0, if the > 30 deg share exceeds X: write node_r/pose_gate_fixposes.txt, move stA0 to stA0_realigned, switch the tree to fixposes and redo the root's stA0 from the input poses.
  - Needs inputs with eo/proj (any refine3D output has them). Resumes from the marker.

## Gate roots (recipe17d_micsup_fpg, SIMPLE_TREE_POSEGATE=0.165 + SIMPLE_TREE_SUPPORT=order, no fixposes key, inputs src_in_fp)
- **The gate fired on both roots:** micnull 31.6 %, mic3 30.4 % of rows turned > 30 deg at stA0.
  - stA0 was moved to stA0_realigned, node_r/pose_gate_fixposes.txt written, and the tree switched to fixposes.
  - Poses then stay at 1.41 deg through stage B.
- **Results identical to the fixposes + order runs:** micnull K* 1 (null PASS); mic3 K* 1 (z -0.65; stage-B groups 12434/12582/13060/12537/6195 as before; state NMI 0.002).

## Summary for the next developer
1. **Main failure: the tree's own re-alignment scatters poses on micelle-dominated small proteins** (fixture 29-32 % turned > 30 deg at the root's stA0; real PfCRT 41 %).
   - Every resulting split is aligned vs misaligned particles, with z 15-73.
   - Masks localise it, so the guard keeps it; without a mask the guard merges it by accident.
   - The real PfCRT tree_umask runs' 20+ leaves are very likely this artefact.
2. **The cross-half EM (<= 19.5 A) loses fine protein differences:** from a perfect truth start it relabels 28-53 %.
   - The 8-deg state is never found. 0 vs 16 deg passes only with the order focus and the truth start (z 4.70).
3. **The order support keeps the inner belt,** so a seeded micelle split passes more strongly under it (z 3 -> 12). The guard keeps micelle splits (0.575).
4. **Tension for any band change:** the RNase P closed differences live at 20-40 A, inside the EM band; the PfCRT protein states live finer.
   - A band change has to be data-derived per node and must keep the 14.6-10.2 A re-test shells in classification (dsnr band and shrinkage failed there).
- **Cost:** ~196 core-h this session (fixture ~70; trees ~122; builds and checks ~4).

## Master posterior merge (scratch build recipe17d_micsup_pm = hetero-tree 5709bc083 + origin/master 9c9ec1105 + my switches; no commit)
- **Roots (micnull / mic3 no mask, plus micnull with PROBBETA=1 and with TOPK=10):** root pose error at stA0 is median 7.3-7.6 deg, p90 120 deg, at every setting (default beta 0.3 included).
- The row's mode itself is wrong for about a quarter of particles at this SNR; tempering does not move the poses.
- **Splits:**
  - Default and beta 1: the split is aligned vs misaligned (pose-error AUC 0.59).
  - topk 10: the split follows micelle radius instead (AUC 0.70; difference energy 0.43 in the belt).
  - Every micnull root: guard 0.255-0.26, merged, one leaf, PASS.

## Soft multi-state reconstruction (scratch builds recipe17d_micsup_soft / _soft2; W/soft/<run>/score.txt)
- **Code** (refine3D refine=cont fixposes=yes, `SIMPLE_SOFT_STATES=beta`):
  - each state is scored against the other half's map at the input pose, with llk = -L P/4;
  - w ~ pi exp(beta llk), pi = the previous iteration's mean weights (every worker reads the part files of iteration t-1);
  - sigma2 comes from the weight-averaged residual;
  - all states are gridded with the weights (calc_3Drec wtab);
  - per part: soft_states_iter<t>_part<p>.txt (row, half, w, llk);
  - diffs: soft_states.diff, soft_states_shrink.diff (scratchpad soft_build/).
- **Result at beta 1:** effectively hard, noise-locked, and it fails every null.
  - mic3 bias60: NMI 0.147 -> 0.010 over 8 iterations; the hard cross-half control (beta 100) is identical.
  - 87 % of weights end above 0.99 at accuracy 0.45; ECE 0.55.
  - The maps fake states: fitted angles 3 / 7.6 / 12 deg while members are near-random.
  - micnull grows a confident split from every start (bias60 vector, seeds, random); random starts do not collapse.
  - Median log-evidence gap is 8 -> 14 nats, against about 1 nat for the ideal classifier.
- **The unchanged split test ACCEPTS the soft-derived noise groups:**
  - micnull seeds groups: re-test z 22, guard 1.025, kept;
  - mic3 bias60 groups: EM regroups 52132 / 4676, z 7.9, guard blind (EM groups).
  - **CORRECTED mechanism:** the soft / soft2 scoring was NOT cross-half. With lp= set, refine3D reads the merged map into both halves (read_mask_filter_refvols l_lpset branch), so every particle was scored against a map containing its own image.
  - The labels therefore depend on both halves' data, and the groups' half maps share the locked noise. That is why the split test is fooled.
  - The soft3 build scores against the true halves; its results are below once they land.
- **SIMPLE_SOFT_SHRINK (soft2) logged gain 1.000 everywhere** (the halves were the merged map), so its six arms are void; two were cancelled by ID.
- **SIMPLE_SOFT_SHRINK as designed:** in the scoring references only, each state's difference from the consensus is kept per shell by max(0, its e/o FSC), with no e/o blend. The maps the split test reads are never shrunk.
  - Six arms (bias60 / seeds / random x mic3 / micnull).
  - Caution (coordinator): this is close to the ancestral shrinkage that failed RNase P r131. Before any tree use, run the r131 closed-subtree root replay (stage B with SOFT_SHRINK in the classification only, then the unchanged split test).
- **e/o blend and the lp= merged-reference rule:** see refine3D_eo_lowres_blend_2026-10-29.md.

## Pose accuracy and findability (fixture, 2026-10-30; W/ceiling_*, W/classify/*, W/trees/*_a1s025_*, *_mixB_*, *_halves)
- **Real pose reproducibility** (pose_compare.py: frame fit, validated on the fixture at 1.41 deg; shiftcmp.py):
  - 56,809 set (pfcrt_core 2_refine3D vs 3_refine3D vs 13_solve3D): disagree at chance (median 126-131 deg, isotropic). Its 4.9-5.8 A FSCs come from lp-set merged-reference refinement. The pfcrt_core / tree_umask* / tree_fixpose* maps and trees are NOT interpretable as conformational results.
  - EMPIAR-10330 stack (SIMPLE 9_refine3D_auto vs the deposited star): 51 % agree within 10 deg (median difference 3.1 deg, ~2.2 deg per run) and 34 % beyond 30 deg. Shifts of the agreeing half differ by 1.0 px at 1.035 A (~0.7 A per run).
- **Findability curve** (ceiling AUC, 0 vs 16 deg, mean-micelle model | hard classifier from bias60, start NMI 0.136):
  - true poses: 0.996 | amplifies
  - 0.47-1.41 deg at 0.31 px: 0.99 | amplifies (NMI 0.25)
  - 2.2 deg at 0.29 px: 0.984 | amplifies (0.235)
  - shift error only, 0.92 px: 0.858
  - the old 1.41 deg / 0.92 px point: 0.820 | decays
  - the measured mixture (34 % random): 0.809 | collapses (0.002)
  - Angle error up to ~2 deg is tolerable; shift error ~0.9 px or a random third is not.
- **The tree at GOOD poses splits the null:** see memory tree-good-pose-null-fails-halves-and-seeding-2026-10-30.
  - Merged-ref stage B: micnull z 10.1.
  - SIMPLE_TREE_HALVES (own e/o half maps; worktree diff, build recipe17d_micsup_fph, gate 15/15): micnull z 4.8, still a fail; 1.41 deg passes; mic3 now partly splits by state (NMI 0.034, leaf angles 3.9 / 8.8 deg).
  - Residual: seeds = flex_pca poles estimated on both halves (OPEN design issue for the whole tree).
- **Soft multi-state reconstruction:** a negative result, see memory soft-multistate-fixed-pose-verdict-2026-10-30.

## EMPIAR-10330 at the deposited poses (W = hetero_tree_20261028/empiar10330; user 2026-10-30)
- **Data:** original unsubtracted stack (16,905, box 300, 1.035 A, EMPIAR_downloads/10330/data/03_Final_Refined_Particle_Stack). The subtracted stack (J407/imported/subtracted.mrcs) is not in the download.
- **Import:** via oritab (rot/tilt/psi = e1/e2/e3, origins px, DefocusU/V/1e4, angast; cs 0.001, fraca 0.1, seeded eo). Every field was verified (<= 1e-5); project e10330/e10330.simple.
- **Gold-standard reconstruct3D at the deposited poses** (no lp, independent halves):
  - all: 7.76 / 3.83 A unmasked, 3.41 / 3.14 A in the 6UKJ region;
  - pose-agreeing subset (deposited vs 9_refine3D_auto within 10 deg, 8,601 particles; e10330_agree): 8.62 / 4.03 A.
  - The deposited poses are good.
- **Frame vs EMD-20806** (pfcrt_core/ref 1.3 A / 240): an x/z axis swap, plus a small rotation / scale still being refined (masks2-4).
- **Roots** (fixposes, maxdepth 0, minpop 250 all / 130 agree, PfCRT keys, CONTRAST=maps, no floor):
  - merged-ref reference: all no-mask K*=1 (z 2.29); agree no-mask K*=2 at z 3.61 (the expected artefact).
  - SIMPLE_TREE_HALVES: all no-mask K*=1 (z 2.47); all 3-class K*=1 (z 3.16); agree 3-class K*=1 (z -0.86); agree no-mask K*=2 50/50 at z 4.15, which is below the residual null (~4.8 at matched pose quality) and so consistent with no states.
  - 6UKJ mask (masks5, corrected frame; SIMPLE_TREE_HALVES): all_6ukj K*=2 at z 5.13. It shows no nuisance association (defocus, astigmatism, CTF FOM / max resolution, SIMPLE corr, pose-agreement flag, micrograph NMI all at chance), and both group maps correlate equally with the consensus (TM / Fab / belt). It does NOT reproduce on the agree subset with the same mask (K*=1, z -0.65). The seed-2 replicate is identical (the root is deterministic).
- **Verdict: no states** in EMPIAR-10330 at the deposited poses. The only split sits at the residual-null level, carries no structure and does not reproduce. No maps were delivered (the user's rule: maps only for a split that survives the checks).

## 16_selection (dibyajyoti PfCRT_smpd_1p3, 44,527 picks, 36 of 50 classes; W = hetero_tree_20261028/sel16; user 2026-10-30)
- **Inputs:** picks matched to the deposited 16,905 by micrograph + position: 8,291 (18.6 %), the rest 36,236. References: the 16.9k consensus halves at the deposited poses, resampled to 240 / 1.3 A (refE_240 / refO_240).
- **Protocol 1, greedy + neigh** (greedy: one iteration at lp 8 against one external consensus half; neigh: 5 iterations, filt_mode=uniform, lpstart 8 / lpstop 5 without lp, so the halves stay gold-standard):
  - A vs B within 10 deg: 12.7 % all, 11.0 % unmatched, 19.9 % matched.
  - matched vs deposited: 27.1 % within 10 deg.
  - reconstruction: 8.00 / 7.80 A unmasked, 6.93 / 5.38 A in 6UKJ.
  - FAIL.
- **Protocol 2, refine3D_auto** (twice; ref_pose_init=cc on refE vs refO; default nonuniform keeps the halves; MAXITS 8, ~64 core-h for both):
  - final FSC0.143 7.26 / 7.09 A;
  - A vs B within 10 deg: 43.1 % all, 40.0 % unmatched (bar 30 %), 56.7 % matched; median 3.9 deg inside, ~2.7 deg per run;
  - matched vs deposited (reference without the particle): 39.7 % within 10 deg (bar ~40 %), median 5.8 deg inside;
  - reconstruction at run A's poses: 9.45 / 8.00 A unmasked, 7.43 / 5.03 A in 6UKJ.
  - Passes at the bar.
- **Deposited-matched vs non-deposited, reconstructed apart** (agreeing picks; 4,699 vs 4,699 of 14,477; no alignment; gold-standard halves; splitmaps/):
  - non-deposited is uniformly worse: 6UKJ FSC at 10-8 A 0.743 vs 0.612; consensus correlation lower in every region (TM .921 / .885, Fab .969 / .939, belt .896 / .847); Delta B +24 A^2; rigid fit 0.0 deg for both TM and Fab.
  - The reproducible (cross-half) difference at 20-15 A, as a fraction of each region's own consensus power: TM 6.0 %, Fab 5.9 %, belt 15 %, rest 13 %.
  - Control, matched from run A vs run B poses: TM 2.2 %, Fab 2.6 %, belt 15 %, rest 19 %.
  - So the belt / rest part is not population-specific, and the specific part is equal in the TM and the Fab. Verdict: worse particles; no helix state.
  - Maps at the top of ~/Downloads: pfcrt_sel16_{deposited_matched,non_deposited}_4699ptcls_{lp10,bm50_lp8}.mrc.
- **Trees** (19,176 A-vs-B-agreeing picks at run A's poses, tree_in/s16_agree.simple; minpop 288 = 1.5 %; fixposes + SIMPLE_TREE_HALVES; fph build; maxdepth 0):
  - s16_6ukj: K* 1 (weakest pair z 2.63 vs z_crit 3.38).
  - s16_c3 (3-class support from run A's reconstruction halves, band 14.2-9.5 A, 27.5 % of the envelope): K* 1 (z -0.05).
  - s16_nomask: K* 1 (z -0.35). TREE HALVES was active in stage B for all three.
- **Verdict: no states** in 16_selection either. The extra picks add worse particles, not another conformation. Cost of the 16_selection work: ~115 core-h (two refine3D_auto runs ~64, the greedy+neigh protocol and validations ~30, split maps and roots ~20).
- OPEN ISSUE (found 2026-10-30, CNGA1 prep): hetero_tree's common_keys hands the tree's pgrp to every flex_pca call, so every OLE tree ran flex_pca with nsym=2. Default unchanged; test switch SIMPLE_TREE_FLEX_C1 (cnga1_tree build). See cnga1_tree_roots_2026-10-30.md.
