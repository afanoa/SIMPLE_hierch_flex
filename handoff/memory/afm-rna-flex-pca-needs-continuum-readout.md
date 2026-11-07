---
name: afm-rna-flex-pca-needs-continuum-readout
description: "On the afm_cryoem RNA data flex_pca's discrete k-center readout fails (76% of particles in one state, rest junk); the heterogeneity is real and peripheral, so use trajectory placement with fat kernels"
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T16:14:54.262Z
---

Run `29_flex_pca` (2026-08-13, 932k particles, 24_refine3D poses+vol, preimage_auto ceiling 16,
shrinkage on): state 10 captured **705,165 particles (76%) at frac_in_support 0.911** and is the
only map that looks reasonable (user-confirmed); the other 15 states sit at 4-63k particles with
frac_in_support **0.009-0.42** — near-empty prior-dominated junk. Kernel bandwidth was floored by
the chi2(ncomp) noise floor at 5.19 while median latent distances are 20-67, i.e. the latent cloud
is 4-13x wider than the per-particle noise floor -- there IS structure beyond noise. Latent noise is
much higher than IgG: shrinkage moved |dz|/|z| by 0.307 (IgG: 0.091). Complete-linkage merge
correctly merged nothing (max disattenuated ratio 0.68 -- the junk maps are noise-distinct).

**Why:** the molecule is a stable core + heterogeneous PERIPHERY (known from the NU local-res
filter), so the latent is a continuum: one dense mode + smooth spread + junk tails. 16 diffusion
k-centers with a noise-floor bandwidth is a discrete-state lens on continuum data -- one center
swallows the mode, the rest land in tails. Same binary works on synthetic (IgG 16->15,
frac_in_support 0.53-0.99; Ribosembly recipe ARI 0.947), so this is data geometry, not a code bug.

**How to apply (next runs, ~45 min each on the dense project):**
- Continuum readout: `npreimages=6 state_axis=1 min_neff=50000 nbins=5 lp=12`, preimage_auto OFF.
  min_neff raises the nn bandwidth floor (fixes frac_in_support), state_axis=1 orders states along
  the dominant component (motion frames), nbins=5 CV-selects bandwidth. All keys are registered.
- Check flex_pca_pc*.mrc for coherent peripheral difference density where NU says heterogeneity
  lives; check state 8 (63k, frac 0.42) as the second-most-real map.
- If eigenvolumes look junk-driven: two-stage -- subset to state 10's 705k (assignments were
  written back to the project) and rerun within the mode.
- 28_flex_pca used the WRONG consensus (23's vol per user) -- ignore its outputs.
preimage_auto + merge is for discrete-state discovery (Ribosembly-like), wrong tool here.
See [[afm-cryoem-1m-run-is-io-bound]] (RESOLVED: dense project), [[flex-pca-thread-over-particles-fix]].
