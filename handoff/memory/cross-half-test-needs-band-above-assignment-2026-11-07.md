---
name: cross-half-test-needs-band-above-assignment-2026-11-07
description: "KEY 2026-11-07 (hier-ppca): a cross-half map test of hard-assigned groups is valid only above every assignment lowpass; in-band, one-conformer null groups read z 53-56 even with swap scoring + view matching; band split (assign lp 15, test 14.6-5.6 A) -> null 1 leaf, qtoy2 classes kept; data fission also passes on the full band; production state_split tests inside stage B's 10 A band (guards carry the nulls); continuum-vs-discrete still undecided"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b40276-4741-4f1d-bb2d-b403a5320b95
  modified: 2026-11-07T19:39:43.623Z
---

hier-ppca side session (handoff `hetero_tree_20261028/handoff/hier-ppca/STATUS.md`, runs `math_reworks/latent/hppca/`).
- qtoynull, lever 4's refined groups: the 5 groups are mutually distinct at z 53-56 (x .30-.33, near equidistant) in 39-10.2 A,
  with swap (other-half) scoring and 300 view x defocus cells. Truth shows no view, defocus, snrq, B or pose-error association.
  Cause: each half was selected by references carrying the other half's noise. Swap scoring alone made it worse (x .35 vs .08).
- Band split: assign (flex seeds and state_eval) below lp L and test above it. qtoynull: 5 -> 1 leaf at L 10 (9.4-5.6 A) and at
  L 15 (14.6-5.6 A, 6 shuffles, max z 1.97).
- qtoy2 at L 15: 8 leaves, NMI .697, closed .916, holo .818 (tree mm2 .780: the PPCA node itself does not beat the tree).
- qtoyflex at L 15: holo .689 (tree .62); the arm is cut in 3; closed lost (.22; tree .39-.47). Nested s1 R2: root latent .187,
  leaf .289, leaf + local .314.
- Data fission is the other valid route and keeps the low band (f2r session, `handoff/fission/`, `math_reworks/latent/fission/`):
  X1 = X + tau Z labels, X2 = X - Z/tau tests. qtoynull tau 1: over-split K 4 -> 1 on lever 4's full band (max pair z 0.46).
- Production (hetero-tree-wip source, checked by the coordinator): `state_split` runs its EM below EM_LIM 20 A and tests from
  kb + 4 shells (focus mask, as the tree calls it; 2 without) to lp 10 A, starting from stage-B labels assigned at LPB 10 A, so its test band sits inside the band that chose
  the labels. The gate nulls pass through `locality_guard` + `effect_floor`; the IgM nulls (10 % 50/50, 5 % 64/36) do not.
- Traps:
  - Per-shell amplitude matching is degenerate for no-signal (junk) maps.
  - Groups under ~2k particles carry no clean-band signal; flag them untestable (half-map corr < .05).
  - 2 shuffles per pair is too few for single-pair nodes.
- Continuum vs discrete is still undecided at map level: bin series show a sorting ramp in-band and have no power above it;
  single linkage on the equivalence graph bridges classes through small mixed groups. IgM 10 % arm info is entirely below 15 A,
  so a band split cannot test it.

**Why:** the user's hierarchical-PPCA model decides states by hard assignment plus a cross-half map test; that test was silently
anti-conservative wherever it shared a band with the assignment.
**How to apply:** any map-level split or merge test (tree split test, merge-up, lever 4) must use a band above every assignment and
pose lowpass, or independent half pipelines (half-specific latents, groups matched across halves). Related:
[[split-test-accepts-quality-split-on-null-2026-10-27]], [[stageb-oracle-references-limit-2026-11-07]], [[discrete-vs-continuous-site-gap-2026-11-05]].
