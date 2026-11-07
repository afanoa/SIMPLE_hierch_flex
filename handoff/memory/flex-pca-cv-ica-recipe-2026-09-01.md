---
name: flex-pca-cv-ica-recipe-2026-09-01
description: "The validated flex_pca recipe (CV+ICA rank, 4-6 iter budget, ridge at low SNR, marching dead on easy data) and the day's measured verdicts on 10028/cryoDRGN"
metadata: 
  node_type: memory
  type: project
  originSessionId: fb1d9e33-0137-47d8-8d28-e3a32330f274
  modified: 2026-09-02T06:24:33.205Z
---

All measured 2026-09-01 on EMPIAR-10028 with cryoDRGN consensus poses
(fwl-c142339, results in /mnt/beegfs/elmlund/afan/10028/10028_simple/flex_h0901).

**The production recipe (halves-first, SIMPLE eo convention):**
paired mod-4 half fits in ONE job (SIMPLE_COV_PAIRED=1, +MERGE=1) at generous
neigs (ceiling), SHORT budget (n_probe_iters=6 or XFSC_STOP patience-3), then
merge -> 2 FULL-SET POLISH iterations (SIMPLE_COV_MERGE_POLISH=2, implemented
2026-09-01, shared-memory only v1; user-proposed, decisively validated: 9-axis
ICA stability 0.847/min 0.184 -> 0.938/min 0.836, ALL 9 axes stable) -> all-N
embedding -> offline CV (principal angles = span) + per-axis cross-half ICA
(latent_ica.py) = k and the physical axes. Validation decisions (rank, axes,
convergence) are frozen BEFORE polish -- polish iterations are unmonitored.
Full-N single fit is retired to a T5-diagnostic role (m1 ships).
Full production run ~11 min at 105k particles.

**Verdicts:**
- Rank on 10028: NO discrete rank. neigs=10 was ceiling-bound (user caught it);
  at neigs=50 the span decays SMOOTHLY 0.97->0.42 with ~30 dirs >0.7 and
  d_signal=21. Core of ~9-13 ICA-certifiable independent axes (6-9 at >=0.85)
  on a continuous tail of small replicable motions (transport picture).
  Deliver at neigs~15 with ICA rotation + stability-ranked axis report; the
  tail stays as regularized variance, not interpreted. The old K*=3 / neigs=4
  calls were artifacts of varimax axis-instability + nested span scoring +
  over-iteration. ICA replaces varimax as delivery rotation. Caveats: polish
  helps marginal axes at small K but at K=50 it spreads variance into MORE
  embedding axes (27->36 healthy) without lifting the 0.7-certified count;
  the in-engine matched-pair S(t) breaks down at K~48 (use span/ICA-based
  statistics at generous K).
- Convergence (overnight trajectory study 2026-09-02, archived bases + offline
  init-DEFLATED principal-angle oracle, 6 arms K={5,10,20,50}x{10k,fullN}):
  oracle t* = 1 in EVERY regime — the data-free init + ONE M-step captures the
  replicable subspace; everything after decays (cost of 6 iters at K=20: 0.94
  ->0.77; K=50: ->0.68; decay steeper at high K and low SNR — the OPPOSITE of
  "more K needs more iters"). Budget: 2-3 (confirmation A/B recipe_it2/3 run);
  40-iter era additionally CREATED dead axes. The old in-engine matched-pair
  S(t) is init-contaminated early + breaks at high K — replace with span-based.
  The period-2 basis oscillation is PROVEN to be the unbalanced cyclic ONLINE
  window (present iff halves>batch: all fullN arms +0.006..0.012 lag2>lag1;
  ABSENT in 10k single-window arms, -0.01). Matches EM theory (O(log) iters) +
  semi-convergence early stopping.
- Cross-fit SSNR ridge (XFSC_REG=1): +0.10 Rule A at 10k, replicated in two
  regimes (box 64 lp16, box 96 lp10). Default ON for small effective N; opt-in
  at full N. Min-blend stays refuted.
- Frequency marching: dead on well-behaved data, ALL flavors — clock neutral
  (identical to 4 decimals), convergence-gated (new SIMPLE_COV_XFSC_MARCH=2,
  implemented+built) WORSE alone (0.76 vs 0.84), invisible with ridge. Cause:
  deepest-crossing finds full band right after warm-up. Keep mode 2 default-off
  for hard datasets (membrane/small/junk-rich) where the band rule would bind;
  battery script reusable (~/run_march_and_k9.sh on the box).
- Distributed paired dispatch (part v5, implemented today on h0901-int2 +
  local mirror): exact at nparts=2, rerun-floor-equivalent at scale, and the
  per-shard ONLINE stripes IMPROVE tail reproducibility (mini-batch effect).
- Box note: fwl-c142339 = 40 physical cores x 2 SMT; budget concurrent jobs
  to ~40-48 threads total.

**Open:** ICA continuum-axis instability is expected (ICA blind to Gaussian
axes); biological readout awaits eyes on flex_h0901/ica_axes/ traversal
volumes (does an axis show 40S rotation?). 10180 transfer pending.

**ONLINE-window limit cycle (2026-09-02):** the fixed-width cyclic window walk
(40000 vs 12623 remainder on 52.6k halves) makes the basis oscillate between
two data-dependent attractors — measured lag-2 self-agreement +0.009 over
lag-1 in BOTH fits. SIMPLE_COV_BALANCED_WIN=1 (implemented, DEFAULT OFF)
equalizes window sizes. USER SEQUENCING DIRECTIVE: hold off on sampling-scheme
changes until the iteration-count question is nailed (overnight trajectory
study); then frequency marching + abinitio3D-style importance sampling
(update_frac random sampling + trailing accumulators — subsumes shuffling and
tail averaging) are the levers to revisit ON CHALLENGING DATA, where they can
actually differentiate — 10028 is too well-behaved to rank them.

**Open diagnosis (2026-09-02): short-budget amplitude anomaly.** At it=2-3 the
delivered latents are ~100x under-variance, but this is NOT prior shrinkage
(prior was weak, 2.6e3) — candidate: early-iteration FSC-Wiener attenuation of
the basis volumes propagating into the projections. The naive Gamma "recal"
(EMBED_RECAL) used the cached MAP SAMPLING precision (Q=A G+ A) as if it were
the posterior precision and COLLAPSED Gamma to 1.6e-2 — default OFF until the
solve conventions are read properly. ALSO: at it=2 the init-DEFLATED merge
match keeps ALL 20 axes >= 0.75 — the 6-iter run's 20->13 cut was partly
pruning axes over-iteration itself had damaged; on this data the prune
threshold is a policy choice (smooth spectrum), not a cliff.

**STATE_ICA verdict REVERSED (2026-09-02, user caught it):** ICA-frame kernel
state placement is a REGRESSION as implemented — the kernel bandwidth is a
chi2-median floor (chi2_median(k), Wilson-Hilferty) that presumes a specific
metric scale; after the frame rotation it produced everything-kernels
(bandwidth 6.219 IDENTICAL across frames = not data-derived; neff 88-93k of
105k per state; frac_in_support 0.95 everywhere) and washed out the
compositional 2nd state and ratcheting entirely. The "balanced occupancies"
were consensus-averaging in disguise. Standing config for STATE deliverables:
default GMM path, raw frame (production20-style — compositional state captured
nicely, per user); ratcheting readout via ICA axis TRAVERSALS (ica_axes
volumes), not kernel states. STATE_ICA stays experimental-off until the
bandwidth rule is derived IN the working frame. General lesson: any
frame/metric change must re-derive every scale-bearing constant downstream
(bandwidths, floors, thresholds) — grep for chi2/floor/bandwidth before
rotating coordinates.

**Day-2 close (2026-09-02):** distributed polish validated (0.993-0.999 vs
shared twin); FSC-doctrine axis weighting live with logged table (at it=2 all
20 axes c=0.79-0.96, w=0.88-0.98 — no cliff, gentle differential shrinkage);
cluster-safe env transport (SIMPLE_ENV_HOOK.sh sourced by all generated worker
scripts; flex master writes it — schedulers don't reliably export env);
regression check PASSED (verify_gmm_br360 = bit-identical occupancy vector
1630/15331/5362/79650/780/2494 to the original full-box run through all code
changes; movement states intact, 780-ptcl state honestly filtered at 6.4 A).
PROFILE (corrected sampler, master+workers separately): master-only phases own
81% of recipe wall at 1.5 busy threads; worker rounds 19% at ~39/80 — the
master, not qsys overhead, is the wall-clock owner of the short-budget recipe.
Fix implemented+built: MASTER THREAD BOOST (master runs its stages at the
workers' combined budget, capped 64, SIMPLE_COV_MASTER_NTHR overrides; phases
never overlap so no oversubscription) — measure on the next profiled run.
Polish round-1 distributed crash fixed (merged eigenvolumes seeded into the
polished namespace before the loop).
