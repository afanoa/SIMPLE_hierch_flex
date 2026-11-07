---
name: coupled-2d3d-reconstructed-tree
description: "~/SIMPLE_coupled, branch coupled-2d3d-rnasep (9eaef62df on master 76c37002b): the RAG coupled-2d3d code rebuilt from the 10-10 + 10-13 handoffs; only subcls_alpha had to be re-created; statically checked, not built yet (2026-10-13)"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T14:33:38.006Z
---

The original coupled-2d3d branch (a72eb3ea1 + uncommitted) exists only on the RAG machine (/home/afan, not batch2,
not the fwl boxes). Rebuilt 2026-10-13 as `~/SIMPLE_coupled`, branch `coupled-2d3d-rnasep`, one local commit
9eaef62df, NOT pushed anywhere:
- base master 76c37002b: the 10-13 patch's known pre-images (simple_fileio, probe_fit_update, state_service) equal it
  exactly, and simple_matcher_ptcl_io = solvent 425241bd1 + patch exactly;
- denoise commander + flex README from 22891675d (base commander imports FLEX_CLS_LABELS_MODE, gone in the copy);
- the nine Fortran files of the 10-13 handoff `code/` overlaid whole (tied module is new);
- re-created: `subcls_alpha` (real, default 100; <0 hard tie, >=1e6 no prior) in parameters/parse/UI, and
  'subclass' allowed in state_placement validation;
- scripts under scripts/cls_expansion_proto/{coupled, coupled/rag10049, subclass_pool, fortran_toy}.

**How to apply:** build on batch2 needs a fresh build dir (new .f90 => cmake), and the args generator trap
([[simple-args-generator-lib64-trap]]) for subcls_alpha. Not on current master: master's 10-12/13 cls_expansion
(view balancing, in-plane off) and particle-io changes would be dropped by the overlay, so porting to master is a
separate step. Recipe and switches: [[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]], the 10-13
handoff note doc/implementation_notes/coupled_2d3d_handoff_2026_10_13.md in that tree.
