---
name: final-estep-ole-real-data-2026-10-24
description: "\"4c\" final K-state E-step on OLE 12707 (both seeds) moves 22-23 % but does not raise cross-seed agreement below the root and erodes the root split (0.986 -> 0.901); FSC flat"
metadata:
  node_type: memory
  type: project
  originSessionId: 5f12aee4-62b6-4518-9284-7a4df90b8e47
  modified: 2026-10-25T03:12:34.197Z
---

2026-10-24, coordinator task "4c on real data". The CAHRA recipe ([[cahra-13991-pose-entanglement-run-2026-10-23]]) was run on OLE 12707 ([[ole-12707-setup-2026-10-23]]):
- Seed 1 = 6_hetero_tree (31 leaves, 14 units). Seed 2 = 7_hetero_tree (26 leaves, including a 35-particle junk leaf; 11 units).
- Scripts in ole_12707/scripts: ole_estep.sbatch, ole_estep2.sbatch (waits for the merge-up), ole_estep_units.py, ole_estep_score.py, ole_estep_2seed.py.
- Output in estep_s1/ and estep_s2/.
- Recipe:
  - start from the root stage-B project with the leaves as states, on m8c55t with the local queue;
  - refine3D runs at box_crop 88 and only gives box-88 maps, so before/after FSC comes from separate reconstruct3D runs at box_crop 200 with postprocess=no;
  - the merge-up units are reconstructed too.
- Cost: 25.4 and 21.3 core-h per seed (80 cores, 16-19 min, 186-210 GB peak).

Results:
- The E-step moves 21.7 % / 23.3 % of placed particles. Sizes even out: the giant leaves shed into small ones (seed 2's r111_s2 goes from 28.4k to 21.1k).
- Per-state half-map FSC is flat: most states change by at most 0.3 A. A few FSC 0.143 crossings jump 1-2 A either way when the curve sits near the threshold.
- Cross-seed agreement, tree -> E-step:
  - root split 0.986 -> 0.901;
  - depth-2 splits at chance before and after (0.55 / 0.52 -> 0.56 / 0.56, chance about 0.53-0.58);
  - units matched 0.308 -> 0.302, NMI 0.341 -> 0.238;
  - consensus about 30 % either way.
- Each seed moves 6-8 % of particles across its own root split in both directions, independently of the other seed. That is what erodes the root.

**Why:** the E-step helped CAHRA because there was one real binary state and an over-broad open leaf. Below a reproducible split, a free K-state E-step over non-reproducible leaves makes near-tie moves; that matches CAHRA's into-open coin flips.

**How to apply:** do not ship an unconstrained final E-step as a general tree stage. If it is kept, constrain it to reproducible structure: no moves across a reproducible split (an E-step per root branch), or E-step states taken from reproducible levels only. Judge it by cross-seed agreement, not by FSC.
