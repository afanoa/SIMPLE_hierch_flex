---
name: tree-verify-components-and-soft-vs-hard-2026-11-06
description: "VERDICTS 2026-11-06 (~300 core-h, softhard/ on batch2): soft/M-est node weights no-go (merged or held-out); ancestral-weight M-est final maps lose to tree-Wiener and dilute closed leaves; refinement mask needed (r131 corner + null z); stA0 / hybrid seeding finds rare states (roots) but oversplits close pairs and loses the r131 r1 corner"
metadata:
  node_type: memory
  type: project
  originSessionId: 9f6eac94-8eb7-4470-9c6e-1b3a2c75c9d3
  modified: 2026-11-06T19:12:10.787Z
---

User 2026-11-06: the tree is "really complicated"; verify every component (knockouts); test M-estimated vs hard refinement with
truncated weights; "do we need mask in the refinement?"; then "maybe M-est needs ancestral borrowing of the weights; per latent weights
too noisy; on the spliceosome M-est improved resolution a lot". Ran autonomously, all of it under batch2
hetero_tree_20261028/softhard/ (scripts also in worktree funny-pike-0046ee scripts/hetero_tree_softhard/, uncommitted).

- **Refinement adds information:** truth NMI seeds -> stage A -> stage B at qtoy2 root 0.35 -> 0.73 -> 0.76, the depth-1 nodes and the detlim
  roots the same way; mic_het root the exception (0.17 -> 0.12). From label files (stage_audit.py). Every stage-B dir keeps the last
  iteration's prob_state_evidence.txt (per-state log-evidence): a free posterior source.
- **Soft node weights NO-GO:** posterior-weighted state member sets are never purer than hard ones on 13 nodes, merged or held-out
  scoring (rare purity 8 %: .296 hard -> .260 soft). Posteriors are overconfident on toys (w .996 right 86-94 %); real 692k stage B
  median max w .89, mass-.99 truncation keeps ~1.9 of 2 states.
- **M-est final maps with ancestral weights:** prototype switch SIMPLE_STATE_WEIGHTS_TABLE in simple_rec3D_strategy.f90 (sparse table
  "row w1..wK", deselects unweighted rows), build softhard/build_mest (fast gate 15/15). qtoy2 vs truth: ancestral/per-particle weights
  near one-hot, no change; dist a=0.4 FSC 6.8 -> 5.0 A but state signal down; tree-Wiener emp2 best on both (4.6 A, dev .694/.149).
  692k RNase P: anc FSC 9.6 -> 9.0 A but closed leaves diluted (r2112 closed-arm .91 -> .81, r1312 1.07 -> .99); emp2 4.6 A, readouts
  intact. The spliceosome gain came because refine3D_states had no other borrowing; the tree already borrows per shell.
- **Refinement mask:** removing it is neutral-to-better on synthetic truth (qtoy2 root .78 vs .76, rare roots purer) but doubles null
  pre-guard z (qtoynull 9 -> 22, olenull 12 -> 20) and r131 r1 loses the closed corner (split merged; stage-B maps still hold 1.02/.47).
  Site mask instead: worse everywhere. KEEP the node mask.
- **stA0-state seeding (no site machinery):** finds the 6.5 % rare state production misses (K* 2, rare .65/.44) and purer rare children
  at 8/10 % (.66/.73 vs .30/.44), nulls pass via the guard, but oversplits qtoy2 (K* 8 root, 4 r1) and loses the r131 r1 corner. A lead,
  not a replacement.

- **Hybrid (user-requested, same day): site seeds + stA0 states as route-A-style extras** (hybrid_seeds.py: n = posterior mass, top-n
  members, skip on FGATE mutual match / floor / half-node / KMAX 8; seed maps via the one-hot importer). Rare roots: 6.5 % found (.72/.34),
  8 % .85/.57, 10 % .85/.66 (prod .95/.30, .95/.44); nulls 1 leaf; qtoy2 root K* 8 with closed class .99/.96 (NMI .74 vs .76) but close
  pairs oversplit (r1 2 -> 4, NMI .29 vs .42; r4 2 -> 3); r131 r1 merged (stage-B 1.03/.50 corner present, pre-guard z 3.8), r2 gives a
  closed-enriched group (purity .111, map .94/.51). Hard gate fails at r1 as forced here; gains were all at roots.

**How to apply:** do not revive soft/M-est node weighting or weight-level borrowing for final maps without a new mechanism; keep pcg_focus
masks; stA0 extras only pay at roots - a root-only hybrid (depth 0, leaving the r131 lineage at depth >= 3 untouched) is the open lead,
needing a tree-code build and the gate battery before any adoption. Related: [[soft-multistate-fixed-pose-verdict-2026-10-30]],
[[latent-first-kmeans-loses-to-tree-root-2026-10-23]], [[klem-integrated-hierarchical-model-vision]], [[rare-state-detection-limit-2026-10-29]],
[[rnasep-closed-state-must-be-retained]], [[spliceosome-mest-refine3dauto-2026-10-02]].
