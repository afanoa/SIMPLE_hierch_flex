---
name: site-split-port-traps-2026-11-07
description: TRAPS met porting the per-site split (xsbg.py / dc_gates.py) into state_mergeup/hetero_tree on 2026-11-06/07; each cost a build-and-gate cycle
metadata:
  type: project
---

Port of the add-only per-site split (`site_split=yes`, commits 043b6bcfc + c4aa7b562 + 378254f33 + 89604e8ac on claude/tree-ports;
worker `simple_site_split.f90`, site path of `cross_em` in `simple_state_mergeup.f90`). Gate 1 = the worker on the 17 saved gate
nodes vs the Python table (handoff/ports/scripts/hetero_tree_mest/xs_gate1*.{sbatch,py,txt} on batch2).

1. Scoring a particle against the site-MASKED group map with the mixing prior and soft weights collapses every site to one group
   in <= 3 iterations. xsbg.py scores the consensus residual against the site-masked deviation of each group from the groups'
   weighted mean, flat prior, hard labels: references = consensus + M (R_g - mu), `l_flat`, `harden`.
2. `cross_em` docks each group map onto the other half's consensus by its own 3-D shift: at a noise-level site the two references
   differ by a consensus-gradient dipole and the EM converges to a translation split both halves reproduce (z 10-35 where the
   Python read <= 3). No docking and no FSC shrinkage on the site path.
3. A node whose site_split worker returned kept step S_SITESPLIT, so `decide` re-ran every poll (147 leaf lines, crash). The step
   returns to S_SPLIT once the summary is read; a worker that left no summary files the leaf; junk-state particles stay out of
   the site children.
4. GAPw's constants (10 A lobe dilation, |B| >= 1330 A^3) only work on the grid they were set on: at the EM's crop (5.9 A/voxel)
   the dilation eats the corridor (B = 0 on the IgM positive control where the Python had 40 voxels). The call runs on the even
   box nearest 5.11 A (maps Fourier-padded).
5. The port's EM has the Python's per-particle evidence at iteration 1 (median |score gap| .16 nats) but it grows twice as fast,
   so it converges within 10 iterations on root nodes where xsbg.py kept flipping half the labels; the root-node factors it
   reports (z 12-38) are reproducible signal and the geometry gates (G2, G3, B) decide.

6. Even on the 5.11 A grid the IgM positive control (site 6) is missed: the split agrees with the Python's at .81 and z1 22, but
   the lobes are blunter and closer (L 14.6 A, B 4 voxels, GAPw .33 vs L 16, B 40, .09). The port's group maps are band-limited
   at the EM's crop (11.9 A); the Python's came from full-resolution reconstructions. Open as of 02:00 EST 2026-11-07: dc_gates.py
   on the port's own maps decides maps vs call.

**How to apply:** when a Fortran EM "finds" far more than the Python prototype, look for coupling between halves (docking,
shared filters) before touching thresholds; run GAPw-type geometry on the prototype's grid. Related: [[hub-hetero-tree-current]].
