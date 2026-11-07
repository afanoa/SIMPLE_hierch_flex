---
name: independent-motions-product-step-2026-10-29
description: "mm3q (3 independent domains, 8 combos): recipe17d tree gives 16 leaves, 1/8 clean combo, f1/f3 re-split in 3 branches, 5 no-motion splits kept; the site layer's 3 root poles combined give NMI .732 and 8/8 combos F1>=.7; proposal = opt-in product step (never RNase P)"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T00:13:35.079Z
---

2026-10-29, independent-motions session (tree_independent_motions_2026_10_29.md in the shared notes; runs hetero_tree_20261028/multimotion;
fixture routes/toy55 simmm3q, make_mm_maps.py / mm_score.py etc.; 115.8 core-h).
- Tree only: 16 leaves (14 genuine), 1 of 8 combinations clean (F1 .89), others .51-.65; f1 and f3 each split in 3 disjoint branches;
  five accepted splits cut one combination in two (not quality, AUC .47-.51) and the merge-up fused none. Best per-motion leaf F1 .28-.54.
- The root site layer already finds one pole per motion (F1 .82-.92, mutual NMI .001). Combining the 3 poles, no pass: combo NMI .732,
  8/8 combos F1 >= .7, per-motion contrast ~2x the best leaf pair. Plus one refine3D pass vs 8 combination maps: NMI .780. A duplicate 4th
  site drops NMI to .700, so an independence filter is required. Leaf-map spatial ICA (+pass) .801 but is region-basis-like: toy only.
- IgM recipe17b tree: the site layer covers 8 of 10 arm angles, leaf ICA 10 of 10.

**Why:** the user (2026-10-29): several independent motions "need to be addressed, perhaps a tree structure is not the best way to represent that".
**How to apply:** proposed design = opt-in key after factorial(): independent root sites -> per-particle product labels
(hetero_tree_product.txt) -> optional one pass -> reconstruct combinations >= minpop. Never run or propose on RNase P
([[no-region-basis-for-rnasep]]). Open: the duplicate-site rule; factor count needs a particle-level held-out test (leaf-map cross-half
rank fails by selection imprint).
- CONTROLS (same day, 2.8 core-h, tool routes/toy55/site_product_readout.py): a naive product puts a fake factor on 7 null runs.
  Gate = pole >= minpop hard members AND decisiveness (hard / soft pole size) >= .5 (set after seeing the controls; held-out check
  pending): real/continuum .81-3.29, fakes <= .40. Gated: nulls 0, micelle 0, rare-state 0, continuum 1, qtoy2 1 (arm swing; product
  NMI .241 vs leaves .706, so no substitute for leaves), IgM 8 of 10 angles, mm3q 3. OLE: 2 motions per seed, both reproduce (|r| .89,
  1.00/.78), cross-seed NMI product .536 vs leaves .32-.33. One re-assignment pass HURTS on OLE (.536 -> .462): use the site poles, no pass.
- HELD-OUT CHECK (offline, gate_heldout.py, 43 new root sites on trees not used to pick .5): 25 real passed, 14 fake rejected, 1 weak real
  rejected (redundant copy), 3 fake passed: 2 weak a3_ph IgM sites absorbed by dependence, 1 real failure = igm_r14 r4_s3 pole 94 % junk
  (decisiveness .71). Fix: drop sites whose pole is mostly amplitude-rule junk (node site_junk.txt; pruned in finished trees). Nulls,
  micelle, detlim: 0 motions. FINAL step: gate (minpop hard, decisiveness >= .5, not junk) -> dependence -> product labels + pole maps, no pass.
