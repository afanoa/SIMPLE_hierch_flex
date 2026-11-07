---
name: flex-pca-basis-rank-over-specified
description: The IgG covariance basis supports an effective rank of ~1 — only component 1 is reproducible across disjoint halves (cos 0.99) and 3-5 are noise; cross-halfset basis cosine is a GT-free signal ranking that works
metadata: 
  node_type: memory
  type: project
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T19:15:59.156Z
---

2026-08-05, neigs=5 IgG run with `heldout=yes`. Cross-halfset principal-angle
cosines between two bases fitted on DISJOINT particle halves, against measured ICC:

| k | cos(theta) | ICC | AUC_gt |
|---|---|---|---|
| 1 | **0.9906** | 0.0276 | 0.4772 |
| 2 | 0.4971 | **0.0296** | 0.5186 |
| 3 | 0.2575 | 0.0002 | 0.4971 |
| 4 | 0.1348 | 0.0000 | 0.4991 |
| 5 | 0.0354 | 0.0012 | 0.4958 |

**The basis is not "unstable" — it is RANK OVER-SPECIFIED.** Component 1 is almost
perfectly reproducible and everything past it collapses. Effective rank the data
supports is ~1, marginally 2; components 3-5 fit noise, as expected when the
covariance being diagonalised is itself 99.7 % noise
([[flex-pca-latent-is-997-percent-noise]]). The two halves' eigenvalue SCALES also
disagree 3x (max 1.88e4 vs 5.88e4), so covariance magnitude is unstable too.

**A GT-free ranking that WORKS.** Spearman(cross-halfset basis cosine, true ICC)
= **+0.600**, vs +0.056 for eigenvalue order and -0.013 for the latent split-half
proxy. This is a DIFFERENT quantity from the split-half rho on the latents that
[[flex-pca-signal-ranking-is-blocked]] shows fails — that one is per-particle
latent correlation, this one is the angle between half-fitted BASES. Caveat: n=5
components, Spearman +0.60 is not individually significant (p~0.28); confirm at
larger neigs.

**But pruning to the reproducible components does NOT help the mixture:** K=5
gives z=+10.82 (top-2), +11.16 (top-3), **+11.94 (all 5)**. Once the per-particle
noise model is calibrated the deconvolution downweights noisy components itself,
so truncation only discards what little they carry. Opposite of the
kernel-regression case, where ordering mattered and truncation still failed.

**HEAD-TO-HEAD vs the linear eigenbasis (run 2026-08-05, same GT contrast, same
cross-halfset estimator, matched rank, Vbar deflated):** flex_pca's 48-component
eigenbasis MATCHES OR BEATS ENM at rank >=30 on all three GT bin pairs
(e.g. b00->b80: eig48(30) +0.200 vs enm(30) +0.150; h00a->h80a: eig48(48) +0.327
vs enm(60) +0.261). ENM wins only at rank 5, on 2 of 3 pairs. **So the
REPRESENTATIONAL claim for the deformation model does NOT hold** — it does not
describe the motion better than a well-fitted linear basis.

**FOLLOW-UP with a PROPERLY HELD-OUT eigenbasis reverses the "ENM wins at rank 5"
part.** basisrep's basis A was fitted on eo=0 only, so scoring it against the eo=1
reconstruction is genuinely held out. At rank 5: eigA +0.256/+0.254/+0.166 and
eigB +0.196/+0.253/+0.130 vs enm(5) +0.130/+0.174/-0.005 on the three GT pairs.
**The held-out linear basis beats ENM at EVERY rank tested — the deformation
model's REPRESENTATIONAL case is fully falsified.** What remains is only the
statistical case (c_basis=0, no rank over-specification), which has never been
shown to translate into better states.

**MOST IMPORTANT: nominally equivalent linear fits differ 5x in explanatory
power** on the same contrast — zhalf5 (neigs=5, 100k) +0.045, eig48 first-5
(neigs=48, 100k) +0.150, basisrep basis A (neigs=5, 50k) **+0.256** — and the best
one saw the LEAST data. So every single-fit basis comparison carries error bars
far larger than the differences being argued over; **average several independent
fits before believing any basis comparison**. That variance-reduction is also the
cheap alternative to switching models at all.

Caveat that cuts the other way: eig48 was fitted on all 100k particles INCLUDING
the ones tested on. Also
`eig5(5)` is markedly WORSE than the first 5 components of the 48-fit
(+0.045 vs +0.150), i.e. refitting at neigs=5 gives a worse basis than truncating a
neigs=48 fit — do not "fix" rank over-specification by simply lowering neigs.
Script: `scripts/flex_pca_basis_headtohead.py`.

**Consequence for the deformation model — STATISTICAL, not representational:** the 38 % basis-error term
([[flex-pca-deconvolved-mixture-works]]) is a SAMPLING defect — it exists only
because the basis is ESTIMATED. An ENM basis is GENERATED from the consensus map
and bead geometry, independent of the particle split, so its cross-halfset angle
is 1.0 by construction and c_basis is identically zero. That removes 38 % of
per-particle latent error and the rank over-specification too (D is chosen, not
discovered) — a quantitative argument for the deformation model independent of the
+0.144 warp-detection result in [[flex-deformation-model-falsified-on-igg]].
