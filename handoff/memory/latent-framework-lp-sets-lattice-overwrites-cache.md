---
name: latent-framework-lp-sets-lattice-overwrites-cache
description: "TRAP: latent_states_v4 sets its lattice from lp (kb = round(box_A / lp), B = 2 kb + 4); setup with another lp REBUILDS and OVERWRITES the shared plane cache (prep_crop has no guard). Install mf_estep's prep_guard; cut bands in code instead"
metadata:
  type: project
---

`math_reworks/latent/absctf/latent_states_v4.py` derives the band limit as `kb = int(round(box_A / lp))` and the lattice as
`B = min(box, max(26, 2 kb + 4))`. `ppinfo_site.prep_crop` returns the cache only if its shape matches `(n, B, B/2)`; otherwise it
rebuilds the planes and saves them to the SAME path. A script that calls `LS.setup` with a different lp (rnasep-misfit tried lp 15.2
to get ring 15: B became 34) would overwrite `tree_latent/cache/r131_lp15.npy`, which other sessions' jobs read. Cancelled in time
on 2026-11-08 (cache untouched).

**Why:** the cache is shared across sessions and takes a full image pass to rebuild; a silent rebuild at another lattice corrupts
every later reader.
**How to apply:** every script that calls `LS.setup` on a shared cache installs the `prep_guard` of `mf_estep.py` (refuse to build or
rebuild); to compare at another band edge, keep lp and cut the band index in code (`mf_nvcmp.py --kbmax`). Flex floors the same band
limit (`projected_model_kfromto`: ring 15 at lp 15), so the framework sits one ring higher than native runs. Related:
[[rnasep-misfit-noise-spectrum-verdict-2026-11-08]].
