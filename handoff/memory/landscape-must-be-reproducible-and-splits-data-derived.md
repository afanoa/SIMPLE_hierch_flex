---
name: landscape-must-be-reproducible-and-splits-data-derived
description: "STANDING (user 2026-10-22): EXTREMA must reproduce (interior labels may be stochastic); the number of splits per node must be data-derived (no K cap / size floor / depth cap), and the conformational landscape must be reproducible AS A WHOLE across independent runs (every level, every map, the continuum's coordinates and occupancy)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T04:59:56.469Z
---

The user, after the r1 audit showed 8-way node partitions disagree across runs (NMI 0.18-0.25): "the number of supported splits needs to
be data derived. the conformational landscape needs to be reproducible as a whole." Earlier: "I don't want to operate with these scripts
because it's slow", "use something that already exists", "we shouldn't even be sampling".

Later the same night (01:45): "I can accept a reproducible conformational landscape as long as extrema are reproducible". Interior
membership may be stochastic. A motion/path-based basin rule is NOT wanted as the general rule ("not everything behaves the same"):
extrema should come from general, reproducible poles (two-tailed soft site rankings along stable directions).
**Why:** a landscape that changes with the seed is not a result; a fixed cap (K <= 8) or a size floor (15k) imposes structure the data do
not support.

**How to apply:**
- No K ladder: the user rejected it ("not very elegant"; on r1 restart stability gave K = 3, but that is no longer the design).
- Extrema come from poles of the stable latent directions: two-tailed soft site rankings, which reproduce across starts at Spearman
  0.99. The interior is sampled.
- Validate any pipeline with at least two independent full runs compared over the whole output (all levels from the root, matched maps,
  continuum coordinates and occupancy), never one run against a reference run.
- Use greedy assignment, and existing SIMPLE modes rather than new strategies.

Related: [[tree-leaves-not-reproducible-across-runs-2026-10-21]], [[overnight-reproducibility-investigation-2026-10-22]].
- 2026-10-28 ~19:30 user: "we should compare conformational landscapes. not individual particles or maps in my opinion" (after
  OLE depth-cut ARI and map best-match tables). Headline metric for run-vs-run = landscape scorer: shared fixed latent space (one root
  flex_pca embedding), states as population-weighted clouds, compare extrema per leading direction, Wasserstein occupancy, span
  coverage; floors = same-build seed pairs + label-shuffle null. ARI / map matching only as diagnostics. Scorer being built by ole-repro.
