---
name: tree-scorecard-and-blur-quality-screen-2026-10-25
description: "score_tree.sh RUN KIND [PAIR] = one-page usefulness scorecard; quality = BLUR of the low-corr group where the map does not change, from UNFILTERED halves; SELF-CALIBRATED cut dB >= dB_ref/2 (node's top vs bottom corr quintile, quality_ref.sbatch) at the SAME BOX (box 88 doubles dB) and SAME FRAME (corr + poses from one project; recB0 frame breaks the 1/2 cut); fixed 55 retired"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-25T20:58:19.564Z
---

2026-10-25, for the coordinator. The user's ask: "a better way to assess if we are doing something useful and if something is worth
doing". Tools are in routes/toy55; pages land in rnasep_klem/scorecard/<run>/.
- `score_tree.sh RUN_DIR KIND [PAIR]` submits score_tree.sbatch on norm, 32 cpus; m33 takes ~10 min. It produces:
  - E: cost from sacct, the jobs whose WorkDir is under the run;
  - A (rnasep): fa_maps10.py --dock on the leaves, plus fa_recall / fa_occ; recall, mutual matches, closed share vs 11.6 %, s1 span
    vs -10.1..+7.5, leaves outside, leaf FSC, corr-AUC leaves;
  - B: split_effect.py on every split's stage-B maps: z, share, local change, ds1 / dclosed via fa, dB;
  - C: validate_tree3 / validate_conf for toys; D: tree_pair_verify for a pair.
- Flags: arm (abs ds1 >= 3 A or abs dclosed >= 0.10); blur; minor (local change < 0.028, the nearest-neighbour deposited median); other.
- QUALITY SCREEN: dB = -4 x the slope of 0.5 log(S_lowcorr / S_other) against s^2.
  - S is the cross-power Re<E_unfil, O_unfil> in the stable region (mask where the 15 A end difference is under 3 sd), 15 A ->
    max(6 A, 2.5 px).
- SELF-CALIBRATED CUT (user rejects a tuned 55): quality when dB >= dB_ref / 2 and >= 3 SE.
  - dB_ref = the node's top vs bottom corr quintile inside dfx deciles, reconstructed as 2 states at the candidates' box
    (`quality_ref.sbatch` TAG PROJ RECSH MASK [CROP]); ~2-3 min at 32 cpus for 692k.
  - Box 160 references: m33 root +177 +- 7, Oct 22 root +148 +- 17, r21 +200 +- 6, r22 (40k) +94 +- 32, cdisc +258 +- 36.
  - dB / dB_ref: quality 0.51-0.64 (site 6 0.64, r3 0.63, s3 0.61, r4 0.51 ON THE CUT); real <= 0.18; root sites 3 / 5 0.29 / 0.28.
- TRAP: box 88 (box_crop, the scorecard's stage-B split maps) DOUBLES dB vs box 160 (root 364 vs 177, r21 443 vs 200). Never
  compare dB across boxes. Below ~40k particles the reference is noise-limited (r22: 3.0 SE at 160, 2.0 SE at 88).
- TRAP: newer builds' node projects (gate replay, base80 stageB/b.simple, recB0/b.simple) hold only the node's rows; `selection`
  refuses them (STK/PTCL2D/PTCL3D counts). quality_ref.sbatch falls back to set_states_binoris.py (patches ptcl3D state slot 23 in a
  copy after checking corr/dfx/e1/e2/state slots against print_project_field on every row).
- TRAP: the quintile reference is FRAME-DEPENDENT (root, box 160): stage B +177/+186, stA0 recB0/b.simple (node_motion's frame,
  identical to motion/cK/m.simple) +298 +- 70, n.simple import -68, stage-B corr with n.simple poses -89. Corr and poses must come
  from one project; dB_ref tracks the quintile corr gap (~55-65 A^2 per 0.01). In recB0's frame site 6 is 0.38: the 1/2 cut (set in
  the stage-B frame) does not transfer; port launch held 2026-10-25 pending a validated cut.
- Option (a) (expected = dB_ref x end-bin corr gap / quintile gap, flag ratio >= 1) FAILED validation (qa_validate.sbatch /
  qa_ratio.py, quality_a/<set>/ratios.txt): flags a strong arm candidate (d80 r1 c4, 1.97), 3/4 unlabelled, 6/7 OLE roots; small corr
  gaps blow the ratio up. Toys (cdisc/IgM) CANNOT calibrate corr rules: quintile references n.s./negative (no simulated quality spread).
  bshare.py (share of end-bin difference = per-shell rescaled mean map) also fails at 1/2. Only labelled quality case = m33 site 6;
  next needs ground truth (toy with per-particle quality spread) or a frame-free calibrator (ptcl2D corr).
- ptcl2D-corr ranking (quality_ref.sbatch RANK=ptcl2D) is NOT frame-stable either: stage-B poses +14 +- 20, recB0 +76 +- 9 for the
  same groups. Quality-screen line STOPPED 2026-10-25; only a ground-truth toy (user's cost/benefit call) could reopen it.
- IgM toy rows (score_tree.sh RUN igm -> arm_score.py: per-angle eta^2 by leaf/unit, span vs truth P10-P90, junk): design
  fulltoy_m77Mdes_igm eta^2 0.233, span 0.60, junk 91 % in 3 junk-majority leaves (r22_s1 97 % junk unnamed), emp2 5.7/10.2 A,
  536 core-h; base fulltoy_m77M_igm 0.221, span 0.72, 79 % in 2, emp2 4.8/6.0 A, 205 core-h -> design does NOT beat base on swing.
- FAILED screens: corr AUC with a shuffled null; cosine to the mean map; stable-region SSNR per particle.
- TRAP: blur from FULL maps is wrong (stage-B full maps are filtered per state: s3 -165 full vs +222 unfiltered).
- m33 row: 27 / 55 recalled, 11 mutual, closed 9.0 %, span -12.0..+7.6 A, 74 splits (arm 21, fixed-55 blur 9, other 45), 28.4 h,
  987 core-h. Against the root's box-88 cut (182), 4 of the 9 blur flags remain (r122, r1122, r122112, r212121).
- Root site 6 of gate/m33_root_m79 is quality, with a 2.1 A arm shift riding on it. Sites 1 / 2 are the arm.
**Why:** usefulness has to be judged against the deposited states and the effect size, and quality has to be told from real states
without the corr-vs-consensus confound or a tuned cut.
**How to apply:**
- Score every new RNase P tree with score_tree.sh against the m33 row.
- Build dB_ref at the box of the maps being judged; never on full maps.
- Tools (all in /mnt/beegfs/elmlund/afan/routes/toy55/; run from rnasep_klem/scorecard):
  - quality_ref.sbatch: TAG PROJ=<node>/recB0/b.simple RECSH=<node>/motion/c1/hetero_step.sh MASK=<node>/mask.mrc. No CROP for node_motion bins.
  - quality_screen3.py: a candidate's dB from bin 1 vs bin 5 unfiltered halves.
  - qa_validate.sbatch: TREE KIND=other|rnasep|cdisc|igm NODES NAME [EXCL]; all-in-one with the corr-gap-scaled ratio_a -> scorecard/quality_a/<NAME>/ratios.txt.
- Related: [[full692-tree-vs-deposited-benchmark-2026-10-22]], [[at-node-switch-swing-call-2026-10-24]].
