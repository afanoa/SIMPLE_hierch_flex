---
name: tree-contrast-and-identity-junk-2026-10-23
description: "hetero_tree contrast + junk (m8c22-m8c32): tree fits its own contrast; amplitude, zero-anchored, identity and smear junk rules all fail on toys or CAHRA; EMPTY ONLY (zero-contrast posterior >= 1/2) is the default; per-particle structural evidence = no-go (cdisc 0.55 nats vs 2.1 needed)"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T19:26:32.440Z
---

2026-10-23 (user decisions relayed by the coordinator; worktree sharp-proskuriakova-b68280, commits not pushed;
snapshots ~/SIMPLE_stree/snap/m8cNN = the build's HEAD).
- TRAP: flex_pca_contrast.txt is 1.0 for every particle (COV_EMBED_CONTRAST_GRID = .false.), so "use flex's contrast"
  changes nothing. The tree fits its own in state_mergeup%fit_contrast:
  - LS amplitude vs the other half's mean map (12 A band) at the particle's own noise (sdn^2);
  - zero atom + log-normal; the posterior-mean contrast feeds the site EM scale and the ptcl3D 'contrast' field (slot 54);
  - that field reaches pftc ctfmats, gridding and PCG (w = a^2, y/a). Unset = bit-identical.
- Junk rules tried, all fail somewhere:
  - zero-anchored: toy recall 0.004 / 0.021;
  - amplitude mixture (>= 2 sd below genuine): takes 9,016 genuine at CAHRA C11;
  - identity test: z 18 at C11 vs toys' 33 / 110; U 0.018 vs 0.014 / 0.033;
  - rotational smear: sigma 0 on toys and C11.
  Toy junk shares the molecule's envelope (map corr 0.89-0.96), so no map-level test tells it from a faint genuine group.
- m8c30 (d09352995): the user's EMPTY ONLY default, junk iff zero-contrast posterior >= 1/2. 0 flagged on toys and CAHRA;
  9,183 at the RNase P root (9,182 of them in its 0.34-amplitude group). The other verdicts are site_junk.txt report columns.
- m8c31 (016ad5f56), step 0: the per-member evidence gap (nats, own map vs matched genuine map, at the member's pose,
  CTF and noise, at the group's mean amplitude) against log((1-pi)/pi):
  - cdisc 0.55 vs 2.10;
  - IgM 6.67 vs 2.00;
  - CAHRA C11 0.61 vs 0.96;
  - RNase P root 8.83 vs 3.30.
  Verdict NO-GO for a universal per-particle junk rule: cdisc recall would be about 0.12 vs the amplitude rule's 0.57.
- m8c32 (94e65051c): fast-tree's node_rows patch, default yes (user decision): step projects cut to the node's rows.
- m8c33 (b6db9a95f): stage A and B one refinement pass each (ITA = ITB = 1; stA0 keeps ITA0 = 3), the user's schedule decision.
- The RNase P root's 0.34 group as maps: stree_port/rnp_root_junkgroup_maps (README). Non-empty 13,122, empty 9,182, random
  genuine 13,122; 52-54 % of the group sits in branch r3 (21 % of random). User: "none of these look like junk" -> NO group
  gate, empty-only stays. TRAP: a pooled map of empties is not evidence of signal (poses aligned to the reference = reference
  bias); the control is phase-randomised images through the same search (not run, coordinator-only).
- Map deliveries match the tree's leaf maps with stree_port/jg_raw.py: the user's states/ maps ARE final_maps/raw
  (tree_wiener 'raw'), and the script reproduces them to 1.8e-7. Set subset states in place with `selection infile=` over
  all master rows (no re-import).
**Why:** the user wants one junk rule that works on synthetic and real data alike; brightness alone flags faint genuine
particles, and toy junk is too close to the molecule for structure to decide particle by particle.
**How to apply:**
- Junk checks replay a node's own site_rank (stree_port/stree_junkval24.sbatch, 32 cores): about 2.5 core-h per toy or
  CAHRA node, 10.6 at the 692k root.
- Score with junkval_score3.py (its column-2 "OLD" now reads the EMPTY flag).
- Steps 1-2 (the pose-marginalised evidence program) need the user's OK.
- USER DECISION 2026-10-23: the junk screen is OFF for RNase P. Every RNase P hetero_tree / site_rank run passes
  `junk_screen=no`. Toys and CAHRA keep the build default (empty-only).

Related: [[factorial-site-layer-port-2026-10-23]], [[margin-weighted-membership-2026-10-23]], [[never-delete-only-read-or-move]],
[[price-new-stages-lead-with-lean-scope]].
