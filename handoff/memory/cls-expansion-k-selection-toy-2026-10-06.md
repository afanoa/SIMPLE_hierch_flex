---
name: cls-expansion-k-selection-toy-2026-10-06
description: "Python toy (scripts/cls_expansion_proto/kselect_toy.py) for choosing cls_expansion's subclass count K and kernel width; held-out mixture likelihood tracks the oracle, Ashman's D stop is broken by construction, the posterior-precision kernel is far too narrow"
metadata:
  node_type: memory
  type: project
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-06T17:46:51.426Z
---

Toy 2026-10-06 (pydev micromamba env; 3600 CTF'd particles/class, tilt-like nuisance, CTF-weighted
rank-8 PPCA, bisection+tied GMM like cls_expansion; scenarios none/discrete3/hinge/mixed; SNR 0.05/0.2/1.0;
lp 12/30 A). Oracle = K minimising ||subclass avg - clean||^2.
- Held-out mixture log-likelihood (Wiener avgs from one half, other half scored with predictive var
  s2 + c^2 var_avg, swapped) picks K within risk x1.00-1.18 of the oracle in every case. Best rule.
- Latent tied-GMM BIC (as logged in simple_flex_pca_gmm) good but overshoots (mixed SNR1: K=6, x1.53);
  ICL undersplits (x1.24-2.02); silhouette biased to K=2 (x1.75 on mixed); even/odd-FRC rule poor (up to x68).
- Ashman's D stop is BROKEN: 2-means halves of a single Gaussian give D~2.7 > 2, so it never stops (K=100
  always). cls_expansion's `sep` readout cannot be used as a stopping rule without a null.
- Oracle K* was 1-3 in all scenarios: at ~3600 ptcls/class, 100 subclasses/class is ~10x worse risk.
- Kernel: at K=50 the held-out likelihood picks h~1.0 (standardised latent; neff ~130-500) and oracle risk
  drops 4-7x vs hard/narrow weights. Real cls_expansion ncls=100 run had neff ~ pop (median 25 vs 20):
  its posterior-precision kernel is effectively hard.
- At lp 30 A in this 48-px/4 A toy the latent carried ~no information (posterior sd > spread): toy scale
  caveat, not a statement about the real 30 A choice.
Next (not done): export 2-3 real parent classes from not/7_cls_expansion and run the held-out K scan.
Related: [[flex-pca-latent-detection-floor]], [[prototype-in-python-before-fortran]], [[flex-pca-must-be-knob-free]].
- BORROWING SIGNAL (same toy, lp 12): rendering each subclass from the model, mu + U zbar_s (zbar_s = mean
  posterior latent of its members), beats the best hard K* average 1.1-5x in risk (SNR 1: 3.0e3 vs 1.56e4
  discrete3, 6.0e3 vs 1.0e4 hinge, 6.3e3 vs 1.46e4 mixed) and stays good at K=100 (hard K=100 is 10-25x worse).
  Per-shell even/odd shrinkage of hard averages toward the model or the parent barely helps (noise estimate
  from ~36-ptcl halves too poor). Model rendering only carries what rank-8 <= lp spans; outside it fall back
  to the parent / shrunk residual.
- BEYOND THE FIT BAND (borrow_highres.py, fit 30 A, data to 8.5 A): extending U to the full band by the PPCA M-step with
  the latents fixed is WORSE than the parent above the band unless shrunk; with per-shell 2F/(1+F) shrinkage (FRC of
  even- vs odd-member extensions) it is never worse than the parent and 10-22% better above the band at SNR 1 (no gain at
  SNR 0.2; FRC weights ~0.00-0.04 in this toy). A full-band fit (higher fit band) helps only where the heterogeneity is
  fine (1 px shift: fit-band risk 6521 vs 7368 at K=100) and is worse above the band. Hard averages lose everywhere at K>=12.
- REAL NOT (2026-10-06, scripts/cls_expansion_proto/not_heldout_bands.py on SIMPLE_FLEXCLS_DUMP dumps of classes 15,
  181, 269; fit on even members, mixture held-out LL of odd members per band, vs parent): heterogeneity IS fine
  (extension FRC weights 0.56-0.72 vs ~0 in the toy). Best everywhere above 30 A at K>=12: latent fitted at 10 A +
  eigenimages extended to full band (fixed-latent M-step, per-shell FRC shrink): 10-6 A gain 16.5-20.6 per ptcl at
  K~100 vs hard 2.6-4.2 (181: -8.1) and 30 A-latent+ext 7.3-11.0. 30 A-only rendering (current Fortran) = parent above
  30 A and visibly 100 parent copies (montage); junk/residual-axis subclasses cannot be rendered. Caveat: Python
  model has no pose-tangent nuisance columns yet, so part of the fine signal may be residual in-plane misalignment.
- PERIPHERY (2026-10-06): NOT region test (not_periphery.py, odd-member reference): extra global rank (16/32) does not help;
  pooled-gated residual helps only at K~12 (core). Toy (periphery_toy.py: rotating core subunit + 4 independent 50%-occupancy
  peripheral sites, CTF-free): REGION-WISE models (core rank 6 + 4 sectors rank 2, placement on concatenated regional
  latents) beat an equal-total-rank global model 2-3x (K=12: core 0.11 vs 0.30, periphery 0.20 vs 0.35 of hard's error;
  K~100: 0.034/0.039 vs 0.055/0.063). Neighbour-borrowed residuals add ~nothing (<=6%). The fix for missing peripheral
  heterogeneity is region-wise latents (user's multi-body idea), not more rank or residual borrowing.
- FOCUSED CLASSIFICATION toy (focused_toy.py, 2026-10-06; dominant core rotation + view/jitter/contrast + weak sites A 50%,
  B 30%, arm): global placement's peripheral purity = chance (0.37-0.41). Per-particle detectability is the bound:
  oracle supervised LDA site A 0.60/0.71/0.87 at SNR 0.05/0.1/0.3 (chance 0.50). At SNR 0.3 focused works: placement on
  data-found hotspot latents ALONE (excess residual variance after the global model, 4-MAD threshold, components)
  purity 0.76 vs global 0.38, arm spread 1.45 -> 1.08 px; two-level global x hotspot only 0.51-0.60. Local tangent
  alignment does not help (and a contrast column inside a peripheral mask removes the occupancy itself). TRAP: the
  even/odd deviation-correlation metric inside hotspots FALLS (0.70 -> 0.34) as peripheral purity rises -> it cannot
  judge peripheral sorting; use assignment stability across halves or held-out likelihood in the hotspot instead.
- FOCUSED COMPOSITE (focused2_toy.py, best version 2026-10-06): hotspots = components of excess variance of the rank-8
  residual (3 MAD) whose residual variance >= 0.5x the variance the global model explains (drops left-over global
  motion; finds 3/3 elements, 0 false at SNR 0.1; at SNR 0.3 the global model absorbs one site -> missed). Per hotspot a
  2-state GMM on a rank-1 focused latent (site accuracy at the supervised ceiling: 0.78/0.74 at SNR 0.1); subclass =
  global state x hotspot states; composite = mu + U8 zbar_g + sum_h mask_h * (members' mean posterior @ per-state
  appearance averaged over ALL particles). Peripheral error vs hard 7-10x lower (SNR 0.1: 0.57 vs 5.57), spread across
  subclasses ~truth (posterior mix; majority-state version overstated it 2.3x). Half-split agreement of hotspot states
  (0.76-0.95) -- RETRACTED as an indicator: see RIBOSEMBLY below, a 1-D threshold split agrees across halves (ARI ~1)
  for ANY direction, state-related or not. Data-found masks did as well as true masks.
- RIBOSEMBLY (~/Downloads/Ribosembly, local, 335k sim. particles, 16 assembly states, star with true poses/CTF;
  ribosembly_focused.py builds a 2-D class from true poses: cone, phase flip, in-plane align, 10 A low-pass):
  all 16 states differ massively -> the global latent sorts them (NMI 0.80 at K=8), focused not needed. Locally
  differing groups ({12,13} 10%, {9,10,11} 8-16%): supervised LDA separates them per particle 91-100%, yet unsupervised
  global/focused/top-30-PC GMM all fail (NMI ~0) except one case (focused 0.33, hotspot 0.54). The state signal sits in a
  single mid-rank PC (5-9). Selecting directions by half-split agreement (trivially ~1) or BIC 2-vs-1 gain both pick
  non-state non-Gaussian directions. Cause: view spread inside the cone (15 deg) creates large non-Gaussian variance;
  at 6 deg the state moves to PC 4-5 and GMM on top-8 PCs reaches NMI 1.00 in one view (0 in the other). View spread
  inside a 2-D class is the main confound for discovering minor heterogeneity -> handle it explicitly or go to 3D.
- VALID HELD-OUT METRIC (2026-10-06): the even/odd deviation correlation is biased whenever the odd particles are
  assigned using the same pixels/frequencies it is scored on (smooth-noise null gave 0.71-0.84 = real 0.75-0.90 on
  PfCRT). Fix: fit and assign on EVEN Fourier shells, score the subclass deviations on ODD shells (independent noise,
  same structure). Positive control Ribosembly (all states, 10 deg cone): true labels 0.86, unsupervised global 0.92,
  random labels -0.01. PfCRT (~/Downloads/03_Final_Refined_Particle_Stack, 16.9k real ptcls, 2-D classes of 333-600
  ptcls at 10 deg, pfcrt_focused.py): every placement ~0, = null -> no detectable heterogeneity in those classes.
- RIBOSEMBLY REALISTIC CLASS (ribosembly_hier.py, 8 deg cone, all 16 states, ~1.7k ptcls, 2026-10-06): two-level
  placement (global tied GMM on 8 PCs inside the particle support, then per global subclass up to 2 focused hotspots
  inside the support, 2-state GMM each) separates locally differing pairs the global placement leaves mixed: view 1
  states 1/2 NMI 0.86 (global 0.00), 3/4 0.54 (0.00), 12/13 0.48 (0.28/0.04), 14/15 0.30 (0.12); view 0 12/13 0.46
  (0.01/0.07), 9/10 0.31 at K_g=16. Overall NMI about equal (0.71-0.81). Some pairs never separate (0/1, 10/11 view 0).
  Shell-split reproducibility of the focused subclasses 0.78-0.86 (true states 0.86, random 0) -> they are real.
