---
name: flex-band-sort-removes-selection-bias
description: "2026-10-03 SIMPLE_COV_EMBED_HP=<A> (batch2 build ~/SIMPLE_bandsort, local worktree ~/SIMPLE_bandsort, uncommitted) sorts flex states on the band above HP only; S4 null bias in <20 A drops from ~580 to ~1 SNR, real axes stay at 9-21"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T00:38:41.470Z
---

**The code (uncommitted):**
- Local worktree `/Users/atifao/SIMPLE_bandsort`, branch flex-band-sort at 90c90a30b (the SIMPLE_region2 commit) plus a patch to `src/main/flex/simple_flex_reconstructor_latent_ops.f90`.
- `projected_model_klo` reads `SIMPLE_COV_EMBED_HP` (Å) once. `highpass_fplane_for_projected_model` zeroes `cmplx_plane`, `ctfsq_plane` and `transfer_plane` for padded indices with h²+k² < (OSMPL_PAD_FAC·klo)². That runs in `prep_imgs4projected_model`, the only producer of E-step and polar planes.
- A band-limited run skips the GPU prep hook.
- State reconstructions (prep_imgs4rec) keep the full band.
- batch2 build: detached worktree `~/SIMPLE_bandsort`, fresh cmake, full build about 4 min. Patch at `~/autostates/bandsort.patch`.

**Test (S4, compose route, n=3 ladders, `rnasep_autostates/S4/hp15`), region-masked end-to-end SNR <20 Å:**

| ladder | full-band sort | 15-8 Å sort (HP=15) |
|---|---|---|
| null directions, region 1 | 526–628 | 0.7–1.6 |
| null directions, region 2 | ~340 | ~0–0.9 |
| real q1 / q2 / q7 | 1179 / 970 / 443 | 21 / 15 / 9.1 |

- The low band becomes unbiased.
- The real coarse change is genuine but far smaller than the full-band sort showed.
- In the sort band (15–8 Å) the null is still biased (49–52 vs real 54–60), as expected.

**Next idea (not built):** cross-band composite state maps. Sort A on band 1 and evaluate band 2, sort B on band 2 and evaluate band 1, then Fourier-splice so every band of every state is unbiased.

Related: [[rnasep-thirds-selection-bias-null]], [[flex-state-maps-selection-bias]].

**Refined thirds (2026-10-03):**
- Setup: band-sorted (HP=15) q2 thirds vs null thirds, refine3D_auto from the consensus start, `~/autostates/as_thirds.sbatch` + `as_r1score_hp.sbatch`.
- Region-1 1−corr at 20 Å: real 0.0042, null 0.0003, random 0.0001. The full-band-sorted equivalents were 0.190 (real) and 0.119 (null).
- At 8 Å, inside the sort band, real 0.154 vs null 0.132, still biased.
- So the unbiased conformational difference is real (14× null) but about 50× smaller than the biased number. Part of the gap is the weaker 15–8 Å sorting (attenuation not yet measured).
- All four refined maps reach 3.50–3.55 Å.
