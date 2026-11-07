---
name: python-latent-needs-ctf-weighted-residual
description: Python prototypes on SIMPLE_FLEXCLS_DUMP dumps - a latent from phase-flipped images picks up defocus (blur/halo axis) and ~1% junk; use c(y - k c P) features + chi2 outlier cut
metadata:
  type: project
---
2026-10-07, NOT classes 15/181/269 (not_render_focused.py): rendered focused averages from a latent of phase-flipped
(sign(c)) images showed a main axis that was a whole-particle blur/halo change (defocus-dependent |CTF| envelope),
not subunit motion; plus a 7-12-particle junk subclass (diffuse blob) at extreme z.
Fix that worked: latent features = c_i (y_i - k_i c_i P) (CTF-weighted residual after the CTF-modulated,
contrast-scaled class mean; k_i also a fitted-not-rendered regressor) and drop |z|^2 > chi2_8(0.999) (~1%) from the
fit and placement. Then the global spread matched the Fortran rendered stack (0.06-0.08).
User correction: a rendered subclass (parent + U zbar) has no population-dependent noise - do not impose minimum
populations; the problem with tiny subclasses was junk, not noise.

**Why:** two wrong stacks before the cause was found.
**How to apply:** any Python latent on the dumps must be CTF-aware; check the main axis of the rendered global
positions in a montage (scripts/cls_expansion_proto/check_render_stack.py) before handing a stack over.
Related: [[cls-expansion-k-selection-toy-2026-10-06]]

Same day, design outcome: the global x hotspot-state GRID (25 positions x 4 states) cost subunit motion and the user
called it inelegant. Accepted form: hotspots = extra LOCAL EIGENIMAGES (rank-1 score of the rank-8 residual inside
the hotspot mask, column restricted to the soft mask) appended to z; 100 subclasses placed in the joint latent,
rendered parent + U zbar (not_render_joint.py). Subunit motion unchanged vs the 8-dim control (90th pct 0.14-0.16;
Fortran stack 0.12-0.13). Fortran port = 1-2 extra basis columns in flex_cls_extend_basis.
