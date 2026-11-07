---
name: batch2-compose-hook-and-worker-names
description: three batch2 traps (2026-10-14): the coupled snapshot build (~/SIMPLE_coupled/snap/9eaef62df) has NO SIMPLE_COV_COMPOSE hook (flex_pca silently runs a plain fit; use ~/SIMPLE_flex2D/build or ~/SIMPLE/build, check with strings); sbatch --export splits comma lists (pass multi-dir sources as separate variables); SIMPLE's distributed workers show up in squeue as "simple_afm_cryoem" under this account, so cancel by working directory, never by name
metadata:
  type: feedback
---

Hit while running the averages-basis compose (`rnasep_coupled/compose_avg`). The env hook is only in builds that
carry `simple_flex_pca_em_compose.f90`: `strings $B/bin/simple_exec | grep -c SIMPLE_COV_COMPOSE` must be > 0
(flex2D build Oct 3: 7 refs; master build Sep 30: 7; coupled snapshot: 0). Without it flex_pca ran a 30-component
paired fit on 692k particles for 38 min before I noticed (no COMPOSE lines, 30 pcs, no flex_pca_compose.txt).
`sbatch --export=ALL,SRC="a,b"` delivers SRC=a; pass SRC1/SRC2 and join inside the script. The master's workers are
separate slurm jobs named simple_afm_cryoem (the qsys default here) with the run directory as %Z; identify them with
`squeue -u $USER -o "%i %j %Z" | grep <rundir>` and cancel by id, since the account is shared.

**How to apply:** before any env-hook experiment, verify the hook string in the binary and look for the hook's log
line within the first minute; cancel wrong runs by run-directory ids. See [[never-signal-simple-processes-by-name]].
