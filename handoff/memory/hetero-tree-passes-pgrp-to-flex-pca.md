---
name: hetero-tree-passes-pgrp-to-flex-pca
description: "TRAP: hetero_tree's common_keys gives every step the tree's pgrp, flex_keys included, so flex_pca runs nsym=N in any symmetric tree (all OLE trees: nsym=2). Test switch SIMPLE_TREE_FLEX_C1=yes (cnga1_tree build, off by default) drops it"
metadata:
  type: project
---

Found 2026-10-30 while preparing the CNGA1 C4 roots (21fac34db). `common_keys` calls `cl%set('pgrp', self%p%pgrp)` and `flex_keys`
calls `common_keys`. So every flex_pca the tree runs (root fflex, SIMPLE_TREE_HALFSEED's per-half flex, the gate refits) gets the
tree's pgrp. flex_pca with pgrp=cN silently sets nsym=N and symmetrises every accumulation ([[cnga1-is-a-continuum-not-discrete-states]]).
- Every OLE tree (pgrp=c2) and every OLE gate has therefore run flex_pca with nsym=2. RNase P, qtoy and IgM fixtures are c1: unaffected.
- Not changed by default, because that would move OLE away from its gates. Test switch SIMPLE_TREE_FLEX_C1=yes drops pgrp from flex_keys
  only (hetero_tree_20261028/cnga1_tree build; diff cnga1_tree_on_halves_v2.diff in the coordinator's shared notes).

**Why:** the user's standing rule is never to pass pgrp to flex_pca; inside the tree it happened without anyone choosing it.
**How to apply:** for any pgrp != c1 tree, decide explicitly whether flex runs symmetric. A default change needs the OLE gates rerun. Open issue for the user (2026-10-30): symmetric-only modes vs C1 discovery,
  to settle before drawing OLE / C4 conclusions.
