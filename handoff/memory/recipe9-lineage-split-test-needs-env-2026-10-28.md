---
name: recipe9-lineage-split-test-needs-env-2026-10-28
description: "TRAP: recipe9-17d production split env is SIMPLE_SPLIT_CONTRAST=maps; RHOMIN=0.03 was dropped 2026-10-29 (guard decides). Pre-guard builds (recipe9-14) still need 0.03 or the null gates split"
metadata:
  type: project
---

2026-10-28: I launched 12 recipe13/recipe14 gates without the recipe9 split env. In those builds the rho floor (rel_effect < 0.03 ->
no split) and the maps contrast are env-only (simple_state_mergeup.f90 SIMPLE_SPLIT_RHOMIN; refine3D SIMPLE_SPLIT_CONTRAST), so the
bare split test ran: qtoynull root split K*=2 (re-test z 14.4 vs z_crit 3.38) = the quality split the floor exists to refuse.

**Why:** the recipe9 lineage keeps the floor/contrast as env switches, not keys or defaults; tree-optim and adaptive-split-merge
export them in their launchers; a gate launcher that does not export them silently tests another recipe.

**How to apply:** every recipe9/12b/13/14 run exports SIMPLE_SPLIT_CONTRAST=maps SIMPLE_SPLIT_RHOMIN=0.03 (and, from recipe14,
SIMPLE_SPLIT_ROUTE unset = choose, =em for the old groups); check run_manifest split_env before reading any null. See
[[check-all-program-keys-before-launch]], [[contrast-switch-kills-real-rnasep-splits-2026-10-27]], [[null-acceptance-majority-leaf-2026-10-27]].
- UPDATE 2026-10-29 (user): the 0.03 floor is DROPPED from production. From recipe16 on (locality guard in the build) launch with
  SIMPLE_SPLIT_CONTRAST=maps and SIMPLE_SPLIT_RHOMIN UNSET (guard on, choose routing). Only pre-guard builds (recipe9-14) still need
  0.03 to keep nulls from splitting. Runs before 2026-10-29 used the floor; check split_env.txt. See [[floor-off-guard-verdict-2026-10-29]].
- UPDATE 2026-11-01: after audit step 1 (claude/tree-dist cleanup) SIMPLE_SPLIT_RHOMIN and SIMPLE_SPLIT_ROUTE no longer exist; only SIMPLE_SPLIT_CONTRAST=maps remains ([[hetero-tree-algorithmic-audit-2026-11-01]]).
