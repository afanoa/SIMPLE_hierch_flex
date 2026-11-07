---
name: tree-user-mask-opt-in-not-default
description: "TABLED design (user 2026-10-27): hetero_tree/flex should take an OPTIONAL user-supplied mask (PfCRT: focus on the helices, not micelle/Fab); never a default solvent mask, it would hurt RNase P / OLE RNA"
metadata:
  node_type: memory
  type: feedback
  originSessionId: b2360bb0-fdc3-4e91-8fc3-27dd729625f7
  modified: 2026-10-27T19:55:36.504Z
---

User, 2026-10-27, after the PfCRT tree showed no step masks out solvent or micelle (focus mask 70 % of the 200 A sphere, 13 % of
it molecule density): "we should probably have an input for a user supplied mask because we need to focus on the helices here. But
this shouldnt be a blanket policy. A solvent mask on the RNAseP or OLE RNA would be detrimental". Tabled for now.

**Why:** on detergent proteins (PfCRT, cnga1) the variance/alignment is micelle- and Fab-dominated and a data-derived Otsu mask
cannot separate micelle from protein; on RNA (RNase P, OLE) a solvent mask clips mobile low-density arms (the closed state).
**How to apply:** when this is picked up, add a user mask key to hetero_tree (and its flex/refine steps), default off; source it
from outside the data (e.g. EMD-20806 docked to the node mean for PfCRT). Never propose a default or automatic solvent mask for
the tree, and never apply one to RNase P or OLE runs. Related: [[pfcrt-tree-run-2026-10-27]],
[[cnga1-solvent-prior-weight-cannot-separate-micelle]], [[rnasep-klem-mask-clips-closed-arm]], [[flex-pca-solvent-mask-dominates-rna-basis]].

**2026-10-29: implemented** as opt-in `tree_mskfile` in recipe17b_umask (see [[pfcrt-core-alignment-2026-10-29]]); default off.
