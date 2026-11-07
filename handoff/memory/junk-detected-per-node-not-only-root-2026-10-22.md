---
name: junk-detected-per-node-not-only-root-2026-10-22
description: "Design (user 2026-10-22 ~13:00, being prototyped by the factorial session): junk is identified at EVERY tree node with one rule, junk PARTICLE = fits no state at its node (per-particle, calibrated null, no fixed threshold); junk STATES (coherent but implausible, e.g. the 107k arm-in-neither-position) are flagged for the user, never auto-removed"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T17:40:58.185Z
---

The user asked "does it make sense to identify junk at only the root?" and agreed it does not, so it is being designed and tested as
independent work in the factorial session.

**Junk particle:** a particle no state at its node explains.
- Junk is only visible relative to a model: the root's coarse reference catches gross junk, and deeper nodes catch subtle junk once
  their states are sharp.
- Test at seeding, where the rest state is formed: score each particle at its pose against every state with other-half references, and
  normalise its best score by what a genuine particle with its CTF / noise should score. The cut comes from a data-derived mixture, not a
  threshold.
- Junk gets label 0 at that node and below.

**Junk state:** a coherent, reproducible but biologically implausible population, such as the user's 107k "arm in neither position"
leaves. Statistics cannot call it junk. Keep it, flag it by quality, and let the user decide.

**Guard:** a rare real state is coherent and passes the halves pole test, so it must never be flagged.

**VALIDATED + DEFAULT ON for new runs (user 2026-10-22 ~15:00; the port implements `junk_screen=yes`).** Spec:
`doc/implementation_notes/per_node_junk_detection_2026_10_22.md` (worktree laughing-carson-4e8ca1); prototype `routes/toy55/junk_score.py`.
- **Score:** best-state amplitude a = NUM/DEN against every converged state, other-half references, 12 A band.
- **Cut:** a 3-component heteroscedastic mixture (extreme deconvolution); junk = components >= 2 sd below the heaviest one.
- **Per-defocus-bin mixture REJECTED:** on real data it flagged 8-34 % genuine particles.
- **simsites:** precision 0.98-0.99, recall 0.65-0.67; the rare 8 % closed state is flagged at the base rate.
- **r1:** no junk.
- **f692st root:** a junk mode of 1.7 % (9,137 particles).
- **The 107k 'neither position':** flagged 2.85 % at the root and 0 at their own nodes, so a junk STATE.
- **Limit:** half-molecule junk is caught about 45 % of the time per node.

**Why:**
- On the simsites toy the new rest state collects essentially all the 20 % junk, so the rest mixes a continuum's interior with junk.
- Junk can also hijack a node's site discovery, which ranks sites by variance share.

**How to apply:** never treat a node's rest state as all interior. Validate any junk rule on simsites truth (rare 8 % closed must not be
flagged) and on r1's rest state.

Related: [[poles-are-minority-vs-rest-k-includes-rest-2026-10-22]], [[f692st-final-deliverables-2026-10-21]].
