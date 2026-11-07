---
name: klem-tree-representation-assessment-2026-10-19
description: "data-derived tree: first-level splits real, deeper nesting at noise level (55-conf); own-half final soft EM is frozen by self-inclusion; cross-half flat re-assignment fixes root mistakes when leaves are clean; recommend tree = proposal engine, output = flat states + similarity graph"
metadata:
  type: project
---

User asked 2026-10-19: "do the children belong to the parents? early mistakes cannot be undone: is a tree the best
representation?" Checks (toy55/: tree_belong.py map test, tree_simmat.py leaf x leaf cross-half similarity of
deviations, tree_flat_estep.py cross-half flat re-assignment with a true-map ceiling; all after the
[[splitbench-rotation-convention-trap]] fix):
- two_50k (2 conformers + junk): every leaf closest to its own branch (3/3); flat cross-half E-step (<=20 A,
  consensus + FSC-weighted deviation) moved 7-9 % across first-level branches, NMI(conformer; leaf) 0.42-0.53 ->
  0.57-0.62 (true-map ceiling 0.75). Equal-N maps over the full band too noisy (NMI 0.16).
- 55_200k lat tree (21 leaves): first level real (39-10 A: 18/21 leaves closer to their own first-level branch, the misses r.4, r.5.1, r.6 = weak clade-2 branches;
  junk / clade 1 / clade 2 blocks), deeper nesting at noise level in the real branches (leaf-vs-sub-branch corr -0.03-0.00; only the junk branch has reproducible sub-structure, corr 0.3-0.4; gt does not record the junk type).
  Final whole-tree pass 21 -> 13 states (6 real + 7 junk), undoing r.5's split, pairing r.1.1.3 with cousin r.1.3,
  merging r.2 with first-level r.4. Three real leaves in three branches (r.1.2, r.2.3, r.5.2; 0-3 % junk) share a
  junk-like >40 A deviation: a shared factor split separately per branch. Flat E-step: 43 % change first-level
  branch but NMI(clade; leaf) 0.22 -> 0.24 only (true clade maps 0.55): leaf maps are clade mixtures.
- splitbench3's final soft EM uses own-half maps: self-inclusion freezes labels (moved 0.1-5 %).

**Why:** answers whether the tree is a valid final description.

**How to apply:** treat the tree as a proposal engine; finish with global consolidation + cross-half flat
re-assignment (replace the own-half soft EM) and report flat states + the similarity graph (simmat), region x state
for independent sites. See [[klem-data-derived-tree-v3-rule-2026-10-19]], [[klem-tree-vs-product-organisation-2026-10-18]].
