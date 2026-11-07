---
name: flex-pca-refactor-branch-state
description: "flex_pca architecture refactor: ONE squashed commit on branch flex-pca-refactor (~/SIMPLE_refactor), rebased onto origin/master 9e0f6e237 as be67a33a4 on 2026-10-01, gate 13/13; ON origin/master as be67a33a4 since 2026-10-01 (fast-forward push); pre-squash history at tag refactor-pre-squash"
metadata:
  node_type: memory
  type: project
  originSessionId: 5112734a-3cc0-4614-86d8-de87e2e517c1
  modified: 2026-09-30T19:27:52.261Z
---

The flex_pca architecture plan (`doc/refactoring_notes/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18.md`)
was implemented 2026-09-30 on branch `flex-pca-refactor` in the worktree `~/SIMPLE_refactor`
(branched from master 6642aaefa, the post-GPU-removal master). 14 commits, phases 0b through 8;
the plan's Section 11 lists commit per phase and the deviations.

**Why:** the user asked for "a nice architecture", implemented on a clean master, with no
builds on the Mac. On the Mac every phase was gated only by a
`gfortran -fsyntax-only` check (`chk.py` in the session scratchpad: stages the flex modules
topologically and syntax-checks the non-flex dependents). Phase 0a (parity fixtures F1-F8)
was skipped because it needs cluster runs, so "behaviour identical" is unverified.

VERIFIED later the same day on batch2: `~/SIMPLE_refactor/build` (fresh configure, gcc 15.2) builds clean,
the fast gate passes 13/13 (job 53502867), and the seeded toy A/B (`/mnt/beegfs/elmlund/afan/refactor_ab/`,
job 53503478, same settings as nogpu_ab/seeded) is byte-identical to master on all 60 outputs; only log
line ORDER differs (per-fit E-step timing/ESTAT/MIX lines now print before the master tail) plus new
per-fit INNER/POLAR timing lines in the paired path.

Phase 9 (2026-10-01, commits f6006c5d4..b3aac7903): simple_flex_pca_em and its 12 submodules are gone.
The fit is `simple_flex_probe_fit` (parent: type + bindings + engine interface) with submodules
`_estep`, `_update`, `_engine`, `_crossfsc` (the polarft_calc precedent: submodules only for ONE
class's methods); the services are plain modules in dependency order util < basis < embed <
probe_fit < pairmerge < fit_driver. 44 flex files. Verified: cmake re-run + build on batch2, fast
gate 13/13 (job 53532329), seeded A/B 60/60 byte-identical (job 53532330). Tooling in the session
scratchpad (p9lib.py parse/assemble, verify.py per-procedure verbatim check, chk.py syntax gate).
Distributed check (2026-10-01, `/mnt/beegfs/elmlund/afan/refactor_ab/distr/`, jobs 53532451/2): master
and refactor masters with nparts=4 nthr=2 on the seeded toy (24 worker sbatch jobs each, ~75 s): all
state maps and eigenvolumes byte-identical, eigenvalues/sig2 agree to 16 digits (last-ulp differences,
embedding 6e-16 relative), i.e. compilation-unit reassociation, not a logic change. Distributed vs
shared-memory differ grossly in BOTH builds (component sign flips, embedding corr 0.43): the known
distribution effect, see [[flex-pca-distribution-changes-results]]. PCG/gridding pairs (2026-10-01, `refactor_ab/pcg/`, jobs 53532565-8, 53532620/1): rec_backend=pcg
maxits_pcg=8, rec_states_backend=pcg, and rec_backend=gridding each byte-identical master vs refactor
(60-61 files). NOTE the toy's DEFAULT M-step is PCG (maxits_pcg=4): a plain run never exercises the
gridding M-step; always add a rec_backend=gridding arm. Every code path the plan names is now covered
except Phase 0a fixtures (never captured).

Phase 10 (2026-10-01, commits 1214d4eff..cdc006607 + docs): library calls for exact duplicates
(gcd, locate_2, euler2m; the rest are NOT drop-in: single-precision pearsn/median, analytic
butterworth, global-state RNG); `simple_flex_pca_records` (flex_selection/flex_fit_model/
flex_latent/flex_state_set) composes the session and 23 wide procedures take the records
(8-24 args -> 2-11); posterior kernels take (fit, ithr) with verbatim bodies in an associate;
obtain_embedding is a 61-line selector (resume / compose / paired), paired axis-weighting +
polish moved into run_flex_pca_paired, the unreachable single-fit branch and
build_covariance_eigenbasis deleted. Tooling: p10lib.py convert() (signature + dummy-drop +
token rename outside strings, keyword-arg aware). Verified again on batch2: build, gate 13/13
(job 53534851), shared/mpcg/spcg/grid byte-identical, distributed last-ulp (jobs 53534882-6).
Lesson: never chain `chk.py | tail && git commit` -- tail's status is 0; test the gate text.

Phase 11 (2026-10-01, f51937ef9): folder layout src/main/flex/{run,fit,states} + root (util, plot,
umap, testers); pure git mv, module names unchanged, the source glob is GLOB_RECURSE so only the
warning-suppression paths in src/CMakeLists.txt moved. Verified: cmake re-run, build, gate 13/13
(job 53536066), seeded A/B byte-identical (job 53536077). The regular reconstruction PCG
(volume/simple_reconstructor_pcg.f90, 3,285 lines) is NOT split either; split that one first if
either is to be split, and let flex copy its shape.

Squash + rebase (2026-10-01): the 14+ phase commits were squashed into one ("architectural refactor of
flex, disentangling numerical and distribution/commander layers"; old history kept at tag
`refactor-pre-squash`) and rebased onto origin/master: first onto c8d888ec1 (23 conflicting files, upstream's
per-file diffs re-applied by hand onto the moved flex files), then cleanly onto 9e0f6e237 -> be67a33a4.
Pushed straight onto origin/master as be67a33a4 on 2026-10-01 (fast-forward, no merge commit); the remote branch was deleted. batch2 recipe: `git push -f batch2 flex-pca-refactor:refs/heads/flex-pca-refactor-in
origin/master:refs/heads/upstream-master-in` (the -f is REQUIRED after a rebase; a rejected push once left
batch2 building the stale pre-rebase tip), worktrees ~/SIMPLE_refactor and ~/SIMPLE_upstream on batch2,
cmake . + make -j16 install, build/fast_gate.sbatch (13/13, job 53539519), A/B harness
/mnt/beegfs/elmlund/afan/refactor_ab/rebase3/{arm,arm_distr}.sbatch with --export ARM,BUILD,EXTRA and
../ab_compare.sh <A> <B>.

**How to apply:** before merging elsewhere, (1) re-run the CMake configure in every build
dir (16 new sources under src/main/flex, model/em_state/rec3D_pcg deleted, em_solve renamed
to posterior), (2) build on batch2 per [[compile-and-test-means-batch2-fast-gate]], (3) run a
short flex_pca job against a master build and compare eigenvalue tables / state maps
before trusting it. The env-knob count in the plan was corrected to 25 real read sites, all
now inside `flex_run_settings%new` (`simple_flex_pca_run_types`). `doc/policies/flex_pca_policy.md`
does not exist. See [[flex-pca-conventions-branch-state]] for the older branch lineage.
