---
name: rho-floor-redundant-with-guard-2026-10-29
description: "TRAP: after the floor drop, launch gates with fixtures/scripts/sst_gate4.sbatch; sst_gate3 still refuses to start without SIMPLE_SPLIT_RHOMIN"
metadata:
  node_type: memory
  type: project
  originSessionId: 3d793c94-9a6d-4a38-bcdd-0083174b8bb5
  modified: 2026-10-29T21:05:32.708Z
---

On 2026-10-29 the user dropped the 0.03 rho floor; the verdict and its numbers are in [[floor-off-guard-verdict-2026-10-29]]. The same day I removed it from my launchers in hetero_tree_20261028.

- **Launchers changed:** f692_r17d.sbatch (the integration launcher), the recipe17b/recipe17e launch_gates*.sh and ole/ole_r17*.sbatch, and recipe14/sst_gate2.sbatch. All of them now export SIMPLE_SPLIT_CONTRAST=maps only and unset RHOMIN, ROUTE and LOCALITY.
  - The f692 and ole launchers refuse to start if RHOMIN is set.
  - Their manifests and split_env.txt record "SIMPLE_SPLIT_RHOMIN unset (rho floor off)".
- **The trap:** fixtures/scripts/sst_gate3.sbatch was left as it was, because other sessions' jobs were running on it.
  - sst_gate3 still requires RHOMIN, so a floor-free launch through it dies at once.
  - The floor-free copy is fixtures/scripts/sst_gate4.sbatch. It also writes rho_floor= into run_manifest.txt and LOCALITY into split_env.txt.
  - r11_child.sbatch now calls sst_gate4.
- **Old versions:** each pre-change launcher is kept beside the new one as *_floor003.*.
- **Left at 0.03 on purpose:** ole_r14.sbatch, igm_replay.sh and sst_f692r9.sbatch, because their builds (recipe14 and recipe9) have no locality guard. The floor_off/ test scripts are also left, set to 0.

**Why:** a new gate launched through sst_gate3 with the new env fails at once, and an old launcher silently re-adds the floor.

**How to apply:** for guard builds (recipe16 and later), use sst_gate4 and the launchers above, and check rho_floor in run_manifest before reading a null. Pre-guard builds keep 0.03 ([[recipe9-lineage-split-test-needs-env-2026-10-28]]).
