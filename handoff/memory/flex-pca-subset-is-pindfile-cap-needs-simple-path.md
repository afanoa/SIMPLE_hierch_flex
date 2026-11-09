---
name: flex-pca-subset-is-pindfile-cap-needs-simple-path
description: "TRAP: flex_pca rowmap= must map EVERY projfile row (it cannot subset; select with pindfile= rows of projfile); sparse_sites cap sub-jobs run $SIMPLE_PATH/bin/simple_exec (default /home/atifao/SIMPLE/build fails 'unknown argument sparse_sites' and the cap runs in-process); coords.txt rows are input-project ids with a rowmap"
metadata:
  node_type: memory
  type: reference
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T17:25:24.615Z
---

Found launching the 50k root sparse pursuits (2026-11-09, `ole_panel/scripts/rootsp2.sbatch`).

**Traps:**
- **`rowmap=` cannot select a subset.** The file must hold one input-project row per row of `projfile`, otherwise "rowmap= must hold one input-project row per project row". To run on a subset, pass the node's full `rows.txt` as `rowmap` and the subset as `pindfile=`: one `projfile` row per line, 1-based. The first launch died at 46 s on all 7 seeds.
- **The cap sub-job needs `SIMPLE_PATH`.** The band ladder's queued cap job runs `$SIMPLE_PATH/bin/simple_exec`. With `SIMPLE_PATH` unset it falls back to `/home/atifao/SIMPLE/build`, which has no `sparse_sites` and exits at once. The main job then runs that cap itself: same results, longer wall.
- **Row ids with a rowmap.** With a rowmap, the `row` column of `sites/coords.txt` holds input-project ids, so `rnasep_closed.py` needs no `--rowmap`.

**How to apply:**
- Export `SIMPLE_PATH=<build>/install` alongside `SBATCH_CONSTRAINT=x6342` in every sparse_sites launch.
- Subset with `pindfile=`.

Related: [[sparse-sites-band-ladder-verdict-2026-11-08]].
