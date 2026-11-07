---
name: flex-pca-lives-as-uncommitted-wip
description: check `git diff --stat src/main/flex/` before trusting flex line numbers; the split landed on origin (2026-08-31) but dirty-tree WIP and divergent local branches remain the norm
metadata: 
  node_type: memory
  type: project
  originSessionId: fd9705ef-2327-4a33-8bc2-088262efeb49
  modified: 2026-08-31T16:04:00.157Z
---

**2026-08-31 status:** `simple_flex_pca_columns.f90` no longer exists — the split
(`simple_flex_pca_em` + 8 submodules) is committed and PUSHED on
origin/flex-pca-split-monolith; see [[flex-pca-em-split-and-syntax-check-recipe]]. The
159-commit flex line sits on local `flex-pca-reproducibility` (unpushed), and the main
tree's only dirty flex file is `simple_flex_reconstructor_latent_ops.f90` (the z-slab
deterministic-insert edit, uncommitted). The habits below still hold — upstream proved it
again on 2026-08-31 when c5dae62b7 removed `init_rec`'s `cached` argument and silently
broke flex compilation on master.

Historical context (pre-split): the live `src/main/flex/simple_flex_pca_columns.f90` was
usually a large **uncommitted** diff in the main tree `/Users/atifao/SIMPLE`, not what
master had. On 2026-08-04 it was +727/-178 against
the last commit, and the differences were behavioural, not cosmetic (packed-CG defaults, an
extra Gram sample buffer). The binaries in `/Users/atifao/SIMPLE/build` and the run logs under
`~/tomotwin100_work/flexpca_run*/` are built from that WIP.

**2026-08-11: the WIP had grown to +2437/-192 across 14 files and local master was 47 behind
origin/master.** Verified with `git log --all -S`, these markers had NEVER been committed on ANY
ref, local or remote — `BATCHED insert_plane_oversamp_multi_scaled` (the batched-OpenMP speedup,
+248 in simple_flex_reconstructor_latent_ops.f90), `STATEREC SPLIT`, `set_merge_default_on` +
COMPLETE-linkage agglomeration, `npreimages=0`/`min_neff=0` discovery defaults, `pcafit`. That is
the 893s->237s speedup work plus auto state-count discovery, existing in exactly one place: the
working tree.

**The SUBSYSTEM is committed** (`a285c845b working flex_pca implementation`,
`654d1f268 Distributed flex_pca ... GMM on manifold`, `3cac8c355 distributed embedding ...`), and so
is `two_gate_state_merge`. Do not say "master has no flex_pca work" — it has the feature; what it
lacks is this delta. The distinction matters when scoping a commit.

**Check both, never one:** `git log --all -S '<distinctive string>'` for whether specific work was
ever committed anywhere, and `git rev-list --left-right --count origin/master...master` for how
stale the local base is.

**How this bites:** a `.claude/worktrees/*` branch is cut from master, so it does NOT contain the
WIP. Line numbers the user quotes come from the main tree and will not match the worktree file,
and the same defect can look different (or absent) in the two versions. A fix landed only on the
branch does not affect the user's runs.

**Do this first** when a flex_pca task quotes line numbers or run-log output:

```
git -C /Users/atifao/SIMPLE diff --stat -- src/main/flex/
```

If it is dirty, read the main-tree file for the real code, and ask before editing there — that
tree holds work in progress that must not be clobbered. A patch built against a copy and checked
with `git apply --check` is the safe way to offer the port.

**Committing one small change is NOT `git add <file>`** — the flex files carry hundreds of lines of
unrelated WIP that would ride along under a misleading message. And `git pull` is usually BLOCKED:
on 2026-08-07 master was 13 behind and 8 of the 12 dirty files were also touched upstream
(`simple_flex_pca_model.f90` among them, with upstream editing the region right next to the WIP), so
git refuses to fast-forward. Upstream edits the flex files too — the WIP is not the only writer.

The move that works, and never touches the dirty tree:

```
git worktree add --detach <scratch>/wt origin/master   # apply the change there, commit, push
git push origin HEAD:master && git worktree remove <scratch>/wt
```

Check first that the target routine is byte-identical at `origin/master` and in the WIP
(`git show origin/master:<path> | sed -n '/subroutine foo/,/end subroutine foo/p'`); if it is, the
change is independent of the WIP and this is clean. The user's local master stays behind afterwards —
say so, since reconciling the WIP against upstream is their call, not a step to do silently.

See [[flex-pca-rename-and-cleanup]], [[tomotwin100-project-setup]], [[no-coauthor-trailer-on-commits]].
