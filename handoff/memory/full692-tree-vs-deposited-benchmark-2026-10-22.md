---
name: full692-tree-vs-deposited-benchmark-2026-10-22
description: "full692 hetero_tree (23 leaves) vs the 55 deposited apo conformers: TRAPS (toy_consensus is a DEPOSITED mean; our maps sit 0.6-3 deg off the deposited frame; the rigid core already differs so absolute CC yardsticks flag every leaf) -> use the paired arm/core test; results"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-23T20:07:53.006Z
---

Analysis 2026-10-22 (batch2 `rnasep_klem/full692_analysis/`, `final10/` and `final10_m/`). Write-up:
`doc/implementation_notes/completed/full692_vs_deposited_ensemble_2026_10_22.md` (untracked, branch flex-falloff-prior worktree).

**TRAPS.**
- `routes/toy55/maps160/toy_consensus.mrc` is the image-weighted mean of the DEPOSITED maps, not our consensus. A "domain
  calibration" against it calibrates nothing.
- The full692 maps (root included) sit 0.6-3.2 deg and 0.4-1.7 A off the frame the deposited conformers were docked into. Rigid-fit
  every map on the rigid core first (fa_maps10.py `--dock`).
- Even the rigid core differs across pipelines: leaf core CC 0.78-0.91 vs about 0.94 between deposited neighbours. So the
  nearest-neighbour CC and the leave-one-out mixture R² put EVERY leaf outside, even mutual best matches.
- Use the paired test: d = arm CC - core CC against the same best conformer. Yardstick: the same d for conformer vs nearest
  neighbour, 5th percentile -0.048 (-0.082 amplitude-matched).
- A core-CC scale factor overcorrects: R² goes above 1, and averages come out "outside".

**Results.**
- **Branches.** r1 closed-leaning 211k; r2 open 185k; r3 open with a weaker arm 151k (lowest core agreement); r4 "neither position"
  145k.
- **Recovered.** 10 mutual best matches, including conf18 <-> r12112 and conf50 <-> r112. 7 of them pass the paired test in both
  variants. Recall is 17/55 conformers (29/55 matched); the open conformers come back best.
- **Outside in both variants.** All of r4, with the deficit only in the arm's swept region (junk-like, as in f692st); and
  r1222, r12111, r12113 (33k, best match conf26 at 0.78).
- **Closed.** 6 leaves at closed-arm >= 0.6, 112k (16 %; 21 % without r4). The deposited set has 7 conformers at that readout
  (11.6 % of particles).
- **Span.** Not wider: outside r4 the leaves span -5.3 to +8.1 A, the deposited set -10.1 to +7.5 A.
- **Stitched embedding** (stitch_final/; union of the 40 node kNN graphs, 61.6M edges): the diffusion map is a three-armed star
  with r3 at the hub; UMAP shows the branches as islands. It largely restates the tree: the only edges between branches come from
  the root latent.

**Swing analysis (user: "are we swinging less than the paper?"; write-up §9; batch2 swing_figs/, landscape/, posetest2/).**
- **Closed end matched at leaf level.** r12112 equals conf18/52 in tip position, reach (97.4 vs 98 A), tip density (2.20 vs 2.14) and
  tip order.
- **Open end stops short.** The open end reaches conf49/47/55 (tip 13-14 A from closed) but not conf43 (26 A). conf43's signature sits
  in the junk-like r4 mixtures.
- CORRECTED (user, 2026-10-23): hetero_tree's refine_keys set prob_athres=0.01 (ATHR). That is the angle threshold for the
  probabilistic SAMPLING, so at 0.01 the sampler is greedy and takes the best candidate. The neighbourhood search around the
  current pose still runs in all three stages (stA0/A at lp 15, box_crop 64; B at lp 10, box_crop 88; 3 iterations each,
  pcg_focus).
- Leaf poses are at a converged local optimum of the whole-molecule objective: a pose test at the tree's settings and the
  prob_athres=10 control both moved 0.03-0.05 deg. They do NOT freeze the angles. I had claimed they did; wrong.
- TRAP: in s47d, automsk=yes promotes envfsc=yes, which rewrites automask3D_state01.mrc every reconstruction. A supplied reference
  mask must be restored before each single-iteration pass.
- **Core-masked alignment (paper-style) at lp 10 on 8-33k particles** scrambles the poses (median 9-18 deg, 90th percentile
  130-143 deg) and erases the arm. Whole-molecule poses are stable (0.03-0.05 deg). So pose compromise does NOT explain a smaller
  swing.
- **No truncation:** mskdiam 240 = 225 exactly; every arm reaches 98 A or less.
- **landscape_sample (m8c5) poles at r1/r2 are not arm paths:** s1 is flat along them. Its held-out end samples are squeezed
  (closed-arm 0.57 at most) because the per-particle coordinate is noisy (fold A vs B 0.59-0.65); the leaves reach further.
**r4 vs conf43 (§9.7):** r4's exclusion is ON HOLD (2D-to-3D session told the user).
- The r4 node matches conf43 in mid-arm position (-27.6 vs -27.7 A), reach and tip, but only as the average of r41 (beyond conf43 at
  half density, best match conf37) and r42 (arm-less).
- No r4 map matches conf43 in arm CC.
**Per-site layer factorial_v3 (2026-10-23, write-up §10; batch2 full692_analysis/factorial/ sites_table.txt, sitecheck.txt):**
- r4: no conf43 pole at normal arm density, and no two-pole conf43-vs-armless site. The five conf43-best poles (12-22k) have masks
  OFF the arm (0.03-0.17 on V), so they carry r4's average arm (mid density 0.35-0.39 vs conf43 0.61); leaves in proportion.
- Root sites: arm r_s1, r_s5; quality r_s6, r41_s6, r41_s7, r2111_s4 (weak: r_s3, r_s7, r2112_s3); candidates r122_s4 (17k,
  extra density beside V) and r1_s2 (33k, density outside the deposited envelope at 76 A).
- HELD-OUT PIXELS (§10.4; factorial/held/, sw_siteheld.sbatch SITE LEVEL = lsamp_score2 + pole re-pick on one pixel set, map
  from the other; 8 core-h for 5 sites): survival r122_s4 0.45, r31_s7 0.42-0.52, r1_s2 / r2111_s3 0.28-0.29, control r121_s2
  0.15-0.29.
  - CORRECTED by §12: the split leaks about 0.31 (fake-template null). So r1_s2 / r2111_s3 are pure imprint.
  - r122_s4 / r31_s7 are partly real (20-30 %). Leak-subtracted they sit inside or at the edge of the conformer range, NOT
    beyond it.
  - Survivors 68 A apart, both next to V. Element NOT named: no RNase P atomic model on batch2 or the Mac.
- CONF43 TEST (§12; c43/, c43null/, sw_c43.sbatch + sw_c43n.sbatch, 9 core-h): UNDECIDED.
  - Held-out conf43 picks (2-20k, arm-only chimeric references) read mid-arm -19..-27 A at 80 % of conf43's density, pass the
    paired test, and are 42-49 % r4.
  - But fake templates reproduce themselves just as well (leak). Conf43's response is only about 25 % above the null.
  - The raw route is blocked: the stacks are 2D-registered, so the poses don't fit raw boxes, and the unregistered stream stacks
    are gone.
  - BAND SPLIT (§12.1; c43band_ctl/, c43band/; sw_bandscore.py + sw_bandheld2.py + sw_bandnull.py + sw_bandamp.py; 5.2 core-h;
    high 14-23 / low 2-10 shells, gap 11-13) is LEAK-FREE: fakes give Spearman 0.02 and about 0 held-out.
  - Result: conf43 is REAL but WEAK. Spearman 0.104; picks +7..+21 SE in the other band (AUC 0.54). Held-out maps lie 8-14 % of
    the way to conf43 (7-19 SD above 20 random draws).
  - That is >= about 2,000 conf43-equivalents, low-band picks 48-54 % r4 / 14-16 % r4111. No clean subset; arm density unresolved.
- ROOT LOW-AMPLITUDE (junk-screen) GROUP on held-out bands (§12.2; lowamp/; 2.9 core-h; maps + README in ~/Downloads rnp_*):
  - The amplitude deficit is imprint: held-out the picks carry 92-98 % of the consensus amplitude.
  - A low-band shape difference survives, view/defocus matched: 2.5-2.7x noise, 38-42 % of its power in the arm region.
  - Picks are r3-enriched (52-65 %) and overlap the port's (a)/(b) only 5-10 %.
  - QUALITY-MATCHED (corr in defocus-decile quintiles, plus view and defocus; sw_lowamp3.py): the arm concentration VANISHES (arm
    share 0.23 = volume) and the residual drops to 1.7x noise.
  - So it is not a separate arm state; it is what low-corr (r3) particles look like.
  - Caveat: corr is measured against the consensus, so it confounds quality with conformational distance.
- Low-corr poles read less closed / more open (closed 0.26 vs 0.34, s1 -4.8 vs -2.5): quality masquerades as an opening arm.
- HOLO (§11; full692_analysis/holo2/, fa_holo2.py): only 8 of 21 holo maps are on batch2 (holo 7-14, EMD-74272-74279, in
  routes/toy55/maps160/holo_docked); the other 13 were never downloaded (~374 MB each). RnpA region P0 = holo-only density
  (9,100 A^3, radius 22-58 A, off the arm): every map of ours reads apo level there (root 0.01; holo 0.45-0.93) -> no protein.
  Holo 8/12/14 have no apo match, and none of our maps matches them better than apo does; r12112/r12111 beat apo by <= 0.013
  raw only. r31_s7 borders the RnpA site; r122_s4 is 24 A from it. Holo maps' core CC 0.76-0.91: use the paired test.
- TRAPS: sites.txt "hard share" = share of the LEVEL's particles with posterior >= 0.5, not pole crispness; pole maps are
  reconstructed from particles selected on their own pixels in the mask, so in-mask excess can be imprint.
**Per-particle quality (2026-10-23, full692_analysis/quality/):**
- r3 is largely a QUALITY split. At matched defocus its root stage-B corr is 0.015 below r1+r2's (AUC 0.12), its pose change
  (dist) is 1.53 vs 0.55-0.61, and it is enriched at low defocus (28 % of the lowest decile vs 18 % of the highest).
- r4 is intermediate (-0.009, AUC 0.27). r42 and r412 are low; r4111 is near-normal, the conf43 candidate.
- The root landscape path t correlates with corr (rho 0.33 in r1), so it is partly a quality coordinate.
- No per-particle ice/CTF score exists in the denoised-stack project; res/res05/cfar there are per-state constants.
Related: [[f692st-final-deliverables-2026-10-21]], [[rnasep-degenhardt-paper-closed-state-and-pipeline]].
