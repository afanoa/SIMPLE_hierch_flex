---
name: averages-basic-sim-placement-failure-2026-10-14
description: "basic sim of flex on class averages (true poses, known mixes): the latent carries the state, the PLACEMENT fails -- standardized all-10-component metric lets heavy-tailed view-outlier components win (GMM peels off view groups); nkern=1-2 fixes pure 100% and paths; spliceosome averages show the same tail signature"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-15T01:10:12.887Z
---

2026-10-14, batch2 /mnt/beegfs/elmlund/afan/toyexp/pv/{avgsim,avgsim2} (scripts avgsim.py make/merge/score,
avgsim.sbatch, avgsim2.sbatch (odd-count padding), avgsim_nkern.sbatch; prototypes in the session scratchpad
avgsimres/place_proto.py, robust_dir.py). The user asked for "basic simulation" after no convincing real-data evidence.

Setup: 2,000 averages simulated directly (ctf=no, SNR 5, uniform random views, TRUE poses) as known mixes w of two
toyc1 volumes: foot (state 5 vs 14, composition) and swing (1 vs 9, region 1 at -20/+20 deg, conformation); w pure
(0/1), uniform U(0,1) or blend Beta(4,4); flex_cavgs=model neigs 10 lp 12, no view centring.

- Pure foot: GMM K2 100 %. Pure swing: GMM K2 51 % (1963/37). The latent separates it linearly at 99.3 % (d' 3.9) with
  1 % view variance: the GMM's best of 13 starts splits off a 37-average group of one view (TV 0.90).
- The state lives in z1 (z1+z2 for blend); the other components are heavy-tailed (excess kurtosis +1 to +50).
  Placement metric is standardized 1/var over all 10, so those components dominate.
- nkern=1 or 2 (state stage only, infile= cache): pure GMM 100 % for both, maps carry 0.96-1.00 of the true
  difference; uniform equal_occ K4 states 0.83/0.66/0.32/0.16 (ideal 0.88/0.62/0.38/0.12), end maps 0.59-0.61
  (all-10: 0.36-0.48).
- Blend (mostly 35-65 % mixes) with all-10: equal_occ ends 0.20 (foot) / 0.03 (swing); untested with nkern.
- Not knob-free yet: signal fraction (deconv) ~1 everywhere, merge cosine high for 7 comps, eigenvalues flat/unsorted;
  t-mixture fixes pure (94-95 %) but not continua; robust t-scatter top PC corr 0.65-0.78 (below z1 alone 0.83-0.89).
  Lightest-tailed components = state components in all six cases.
- Spliceosome averages latent: comps 1-2 lightest tails (+1.3/+1.7), rest +2.5..+12.4 -> nkern=2 real test
  (full320/flex_avg_{equal_occ,gmm}_K5_nkern2) launched.
- Traps: simulate_particles refuses an odd nptcls ("odd number of projections directions"); flex refuses a project
  copied from a flex run dir (carries state labels) -- resume from the original project + infile=embedding.

- Paired sims (avgsim3: 1,000 sub-classes x even/odd draws): every latent component reproduces between the pair
  (reliability 1.00) -> the extra components are deterministic view-specific structure, not noise; GMM K2 fails for the
  swing under every variant (t-mix, trimmed, uniform background, residual-trimmed, reliability-weighted); equal_occ path
  states reach the blend ceiling (ends 0.37 foot / 0.30 swing vs ~0.40 ideal; maps 0.21-0.30).
- Spliceosome real averages (34-class expansion): even/odd agreement of the sub-class latent z1 0.19, z2 0.11 (noise),
  z3 0.44, z6 0.38; regressing out the view offset within the parent (quadratic) leaves z3 0.29, z6 0.35, z1 0.09:
  the reproducible part sits in HEAVY-tailed components, so "lightest tails" (sim rule) picks wrong on real data; the
  averages carry little reproducible heterogeneity whatever the placement (best label kappa 0.53, one 90 % state).
  nkern=2 on real data: end-state difference unchanged (r 0.73-0.42 vs 0.79-0.41 all-10).
- Next tested: view-tight 2D (100/200 classes) -> kernel expansion, gate 1 = particle-flex state enrichment in the
  sub-classes (old x1.2-5), gate 2 = even/odd latent agreement (scripts chain2d.sbatch, enrich2.py, evgate.py).

User: GMM works well on particles; a separate placement mode for expanded averages is acceptable.
Related: [[averages-route-real-data-2026-10-14]], [[averages-flex-pv-toy-validated-2026-10-14]].
