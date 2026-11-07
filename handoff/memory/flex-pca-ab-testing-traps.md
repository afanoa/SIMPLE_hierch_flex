---
name: flex-pca-ab-testing-traps
description: "Four traps that silently invalidate flex_pca A/B comparisons — worker binaries, state labels, nparts=1, and cherry-picked reference volumes"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e31f34db-3e71-4d09-ab25-71f10f48fa4c
  modified: 2026-08-07T13:39:36.346Z
---

Four ways a flex_pca comparison can look conclusive and be worthless. All four bit on 2026-08-06.

**1. The project pins the worker binary.** `sp_project`'s compenv stores `simple_path`, and
`distr_simple_script` launches workers from `$simple_path/bin/simple_private_exec` — NOT from
whatever binary you invoked. So `nparts>1` runs of an old build silently use the CURRENT tree's
workers, and the workers do the SNR, column accumulation, probe EM and embedding. Two "old code"
controls were invalidated this way before it was noticed.
**Why:** the master binary you type is only the orchestrator. Exact mechanism:
`simple_qsys_env.f90` `qsys_env%new` does `if( params%projfile /= '' )` -> read the compenv segment
from the PROJECT and use its `simple_path`; **`$SIMPLE_PATH` is consulted only when there is no
projfile**, so exporting it before the run does nothing.
**How to apply — the clean fix (no data touched, no nparts=1):** copy the project into the run dir
and repoint the copy, then run against the copy:
`cp orig.simple $OUT/ && SIMPLE_PATH=<new build> simple_exec prg=update_project projfile=$OUT/orig.simple`.
`update_project` rewrites only projinfo+compenv (`read/write_non_data_segments`), and IgG/Ribosembly
stack paths are absolute, so the cwd rewrite is harmless.
**Pre-flight, always, before trusting a distributed run:**
`strings <projfile> | grep -o 'simple_path=[^ ]*'` and after launch
`grep -o '^/home[^ ]*simple_private_exec' distr_simple_script_1`. Do this BEFORE the run, not after.
Do NOT reach for `nparts=1` here — that is trap 3 below.
**Worse than that:** `make` writes to `build/production/` and does NOT refresh `build/bin/`, while
`$PATH` resolves to `build/bin/` — so the MASTER goes stale too. On 2026-08-06 `build/bin/*` was
20+ hours behind `build/production/*` and an entire day of source edits had never executed; a
source revert to "the good code" was therefore a no-op. Always
`cp -f build/production/simple_exec build/production/simple_private_exec build/bin/` after a build,
copy BOTH (a master that stamps a new cmdline key kills a stale worker), and `md5sum` to confirm.
Never trust the `**** SIMPLE Git Commit:` banner — it is stamped at CMake configure time.

**2. State indices are arbitrary k-means labels.** Comparing "state 3 to state 3" across runs is
meaningless. Doing so made runs look flat that were not: index-matched pairwise correlation said
0.77-0.80, while the label-invariant MOST-SPREAD-TRIPLE said 0.61 against the reference's 0.57.
**How to apply:** always score over all C(n,3) triples and take the most spread, or match states
optimally first. Never compare by index.

**3. `nparts=1` is not a neutral control.** It made 19 of 20 IgG states collapse to featureless
balls with code that is fine at `nparts=10`. Convenient for isolating a binary (trap 1), but it
changes the answer, so never treat a `nparts=1` result as representative.

**4. Low pairwise correlation does NOT mean "states differ".** A degenerate blob decorrelates from
everything and scores as the best result. Score TWO things: `spread` (mean pair corr of the
most-spread triple, lower = states differ) AND `cmin` (the worst state's fraction of positive mass
in its densest 1% of voxels, higher = still a molecule). On IgG: reference 0.574/0.345, good runs
~0.61/0.39, the blob run 0.100/0.222. Scorer lives in the session scratchpad as `score.py`.

Also: correlation-to-consensus is NOT a usable quality metric across projects — the Aug-5 reference
volumes score only 0.22 against IgG-RL's consensus because they were built against a different
refinement's consensus. See [[flex-pca-speedup-verified]].
