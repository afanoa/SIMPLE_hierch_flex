---
name: flex-refine3dstates-gridding-and-knobs
description: "refine3D_states forces rec_backend=gridding for its flex init (no solvent mask possible); preimage_auto is unreachable through it; min_state_frac overflows at high nstates; n_probe_iters=2 beats the default 5 on a masked basis"
metadata:
  type: reference
---

Verified 2026-09-19/20 on master 34e69f0c4.

**refine3D_states forces gridding for the flex init.**
`simple_commanders_refine3D.f90:1228` in `run_flex_pca()`:
```fortran
! the gridding M-step for the initialization: no PCG kernels and no solvent support needed
call cline_flex%set('rec_backend', 'gridding')
```
So `flex=yes` inside refine3D_states ALWAYS solves the basis on the spherical
`mskdiam` support — `pcg_mskfile` is inert there. Standalone `flex_pca` is different:
its commander defaults `rec_backend=pcg`. To get a solvent-masked basis you must run
`flex_pca` standalone, then `refine3D_states flex=no` on that project (it picks the
labels up via the `NSTATES FROM PROJECT` branch). The user's view: the *refinement* is
what produces good states, so masked-flex -> flex=no refinement is the combination.

**preimage_auto is not reachable through refine3D_states.** It hard-requires
`nstates >= 3` on the cline and passes it down as the flex `npreimages`, so the
`AUTO_NSTATES=32` branch (which needs `npreimages` undefined) can never fire. It DOES
read the delivered count back (`nstates_flex = get_n('state')`), so the wiring gap is
small: skip the `npreimages` set and relax the precondition. Standalone `flex_pca
preimage_auto=yes min_state_frac=0` works today.

**min_state_frac overflows at high nstates.** refine3D_states defaults it to 0.1 under
`flex=yes`; the check is `nstates_req * nint(min_state_frac*nptcls) > nptcls`. At
nstates=10 on 932,496 ptcls: 10 * 93,250 = 932,500 > 932,496 — fails by 4 particles.
Use `min_state_frac=0.05` at 10 states. It is also incompatible with `preimage_auto`
(THROW_HARD) so auto mode needs `min_state_frac=0` and therefore has NO population floor.

**n_probe_iters=2 beats the default 5** on a masked basis (trpm4, 153k): final
`|cos|` 0.942 vs 0.819, same 3 macro-clusters, 16 min vs 25. Consistent with
[[flex-pca-optimal-iters-depends-on-rank]] — optimal iterations fall as the basis sharpens.

**SLURM job state lies.** The generated scripts pipe through `| tee`, so a hard STOP in
simple_exec still reports `COMPLETED 0:0`. Always grep the log for `ERROR!` / `NORMAL STOP`.

**particle_sieving is not deterministic**: identical coarse parameters gave 79,686 vs
116,836 coarse rejections on two runs of the same flipqr data (~5% swing).
