---
name: flex-pca-amplitude-mstep-omission-2026-10-27
description: "flex_pca amplitude: M-step omits a, a^2 (fix alone ~no gain); delivered embedding at a=1 -> fitted-contrast embedding removes quality from the latent (qtoy snrq R2 .17->.01) and raises seed ARI; the sharp-proskuriakova fix (a, a^2) is guarded by l_sp_req so the production plain/tied path is NOT corrected"
metadata:
  type: project
---

2026-10-27, read of worktree sharp-proskuriakova-b68280 src/main/flex/fit/simple_flex_probe_fit_estep.f90 (fit_estep_solve_stats).
- a = clamp(myv/e_mm, 0.1, 5) enters the posterior (a^2 G / sig2) and the residual y - a T mu (subtract_mean_banded).
- The model y = a T (mu + W z) needs a E[z] and a^2 E[zz'] in the basis normal equations. The local fix scales zbatch/dens
  only `if( fit%spec%l_sp_req )` (shared/private arms); plain and tied-MCFA production fits still pass unscaled moments.
- Investigation doc (~/Documents/SIMPLE_Investigation_Changes_and_Test_Suite_2026_10_27.docx): rank-one counterexample loss
  3.91 -> 3.395 corrected vs 6.608 uncorrected; no native acceptance yet. Separate open issue: training 1/sig2 vs evaluator 2/nu.
- I traced only the E-step side; the doc's audit says no M-step consumer (gridding, PCG) supplies the factors.

**Why:** flex_pca is the user's code and feeds every tree's fflex latents; the omission biases W whenever a spreads.
**How to apply:** raise it with the user before any change; never widen the guard unasked. Related: [[flex-pca-deconv-precision-mismatch-flag-2026-10-22]].

A/B RESULT 2026-10-27 (batch2 hetero_tree_regression_20261026/amp_mstep_v1: clean d3513ff70 vs d3513ff70+patch
`.claude/patches/flex_amplitude_mstep_on_8c1b66f31.patch`, identical flags, fast gate 15/15; frozen root fflex replays of
r1211 / qtoy2 / OLE 6_ node_r, seeds 1+2, tree part counts; score_summary.json; 28.1 core-h):
- sd(a) is NOT small: r1211 0.16, qtoy2 0.24, OLE 0.26 (mean a 0.98-0.99).
- The fix moves the delivered W*z field 10-14 % RMS (cos 0.99-0.995), 4-18x the seed spread; labels reshuffle (ARI 0.43-0.68).
- Truth gains are tiny: qtoy closed AUROC 0.976->0.978, closed_arm R2 0.604->0.613; r1211 closed AUROC 0.709->0.715. Small
  costs: qtoy junk AUROC 0.876->0.862, quality (snrq) R2 0.123->0.130, r1211 best closed k-means label F1 0.37->0.32-0.34.
- Verdict: correct algebra, real change, no fit-level win. Not promoted. Hypothesis: delivered embedding uses a=1
  (COV_EMBED_CONTRAST_GRID off), so amplitude leaks into the latent at delivery; next arm = fix + amplitude-aware embedding.
- TRAP: d3513ff70's qsys_env takes simple_path from the RUNTIME SIMPLE_PATH, so the old "project pins worker binary" trap
  no longer applies there; a grep for the build path in the project finds nothing (my pin check false-failed once).

FACTORIAL 2026-10-27 (same harness; arms base / amp = M-step fix / emb = embedding at the probe E-step contrast
a = clamp(<Tmu,y>/<Tmu,Tmu>, 0.1, 5) instead of a=1, `.claude/patches/flex_embed_fitted_contrast_on_d3513ff70.patch` /
ampemb = both; builds amp_mstep_v1/v_emb, v_ampemb; all fast gates 15/15):
- The EMBEDDING fix is the valuable one. qtoy2: quality (snrq) R2 of the latent 0.169 -> 0.013, closed_arm R2 0.604 -> 0.661,
  s1 R2 0.603 -> 0.692, closed AUROC 0.976 -> 0.979, junk AUROC unchanged; cross-seed label ARI up on all three datasets
  (qtoy 0.74 -> 0.95, r1211 0.86 -> 0.96, OLE 0.95 -> 0.97). Costs: k-means best-closed-label F1 qtoy 0.65 -> 0.47,
  r1211 closed AUROC vs the historical reference 0.709 -> 0.697 (reference labels are not independent truth).
- ampemb: same latent gains as emb but k-means labels become seed-unstable (ARI 0.20 r1211, 0.37 OLE, 0.51 qtoy).
- Delivered contrast with emb: sd 0.16 (r1211) / 0.24 (qtoy) / 0.26 (OLE).
- Relevance: quality-driven splits (RNase P r3, m33 root site 6) are a standing complaint; emb removes quality from the latent.
- Next (launched): tree-level base vs emb on r1211 enriched (maxdepth 2), qtoy2 and OLE root+children (maxdepth 1), seeds 1-2,
  qtoy also ampemb; dirs amp_mstep_v1/trees/<ds>_<arm>_s<seed>, launcher tree.sbatch.
- TREE LEVEL (13:05, amp_mstep_v1/tree_score_summary.json; qtoy root+children maxdepth 1, r1211 enriched maxdepth 2):
  r1211 s1: emb best closed leaf F1 0.451 -> 0.574 (P .69/.66, R .33 -> .51), 5 -> 4 leaves (s2 pending).
  qtoy s1+s2 (consistent): closed recall .92 both, precision .35 -> .26 (F1 .50 -> .41), struct NMI .51 both, leaf-quality
  association already tiny in both (NMI .004, MI|conformation .024). Mixed; OLE full trees + nulls decide.
- TRAP: an ssh output filter `grep -v "^  "` (meant for the login banner) also strips indented JSON lines; fetch JSON with
  2>/dev/null and parse from the first '{' instead.
- 14:00 seed 2: r1211 emb closed best-leaf F1 base .45/.38 -> emb .57/.61 (R .33/.31 -> .51/.58, P .69/.47 -> .66/.64):
  consistent gain on real RNase P across seeds; qtoy2 slightly worse (F1 .50 -> .41, both seeds).
- NULL TREES (qtoynull, published build, no contrast fix): BOTH arms split the null root by particle quality
  (base: 33.5k/26.5k, snrq means .85 vs 1.33, AUC .12; emb: 40.9k/19.1k, AUC .86); base then splits r1 (AUC .74) and r2.
  The embedding fix does NOT stop quality splits; they arise in refinement + the split test -> SIMPLE_SPLIT_CONTRAST=yes
  trees (adaptive session) are the test. amp_mstep_v1/null_root_check.json.
