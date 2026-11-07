---
name: relion-multibody-benchmark-values
description: "Nakane et al. 2018 eLife 7:e36861 numbers for 10180 (4 bodies, consensus 4.3 A on 327k, SF3b 5.1 -> 4.4 A) and 10028 (3 bodies, head 4.0 -> 3.7 A, PC1 bimodal) to benchmark SIMPLE region flex"
metadata:
  node_type: memory
  type: reference
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T04:42:28.769Z
---

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6005684/ (the elifesciences.org page does not render for WebFetch).

**10180 spliceosome B-complex:**
- 327,490 polished particles, box 320 at 1.7 Å/px; consensus 4.3 Å.
- 4 bodies: core, foot, helicase, SF3b.
- Round 1 per body: core 3.8, foot 4.0, helicase 4.5, SF3b 5.1 Å.
- Round 2: SF3b 4.4 Å, the others unchanged. Masks were redrawn from the round-1 body maps, and the round-1 body maps low-passed to 7 Å were the starting references.
- Priors: 10° angular, 2 px translation.
- PCA: two components ≈30% of the variance, unimodal (continuous). PC1 = SF3b rocking over the core; PC2 = helicase and SF3b rocking together (towards Bact).

**10028 Pf80S:**
- 105,247 particles; consensus 3.2 Å.
- Bodies: LSU 3.1, SSU 3.2, head 3.7 Å (head local resolution 4.0 → 3.7 Å).
- PCA: two components ≈30%. PC1 is bimodal; classes on it refine to 4.4 and 3.2 Å.

**Not in the paper:** motion amplitudes (degrees or Å) and EMDB codes for the body maps.

**Our 10180 work** uses 46,783 particles (5.91 Å consensus, `10180/consensus46k_refine3D_auto`), so compare improvements relative to our own consensus, not absolute values.

Related: [[relion-multibody-and-symexpand-facts]], [[auto-regions-recover-ribosome-subunits]].
