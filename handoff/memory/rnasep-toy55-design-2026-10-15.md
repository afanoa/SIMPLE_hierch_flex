---
name: rnasep-toy55-design-2026-10-15
description: "RNase P realistic toy (toy10): particles simulated from 10 deposited Degenhardt conformers (EMDB, mirrored into our hand, box 160 @ 1.464 A) with real refined poses + CTF, image counts as weights, bfac/snr calibrated to the real band profile (signed CTF, snr grid 3e-4..1e-2), 20 % junk, refinement from perturbed poses; toy10c adds two compositional variants; scripts, batch2 layout, traps"
metadata:
  type: project
---

User (2026-10-15): "the data supports a far richer conformational landscape yet any attempt to model it ... has
failed ... our tests are not thorough"; then "we just need 10 conformers instead of wasting time"; then the standing
redirection: "I do not care about replicating their workflow. I want a LATENT VARIABLE approach that can model large
motions and conformational and compositional variability" (see [[toy10-latent-night-2026-10-15]]).

Design (`scripts/cls_expansion_proto/toy55/`, batch2 `/mnt/beegfs/elmlund/afan/routes/toy55/`):
- Conformers: `rnasep_conformers10.csv` = conformers 1,10,18,27,28,37,46,48,49,51 (EMD-74131/74201/74209/74230/
  74231/74246/74259/74132/74133/74074; closed = 18, 51 = 14 % of the mass). `process_maps.py`: Fourier-crop x2 to
  1.464 A, 160 box around the density centroid, z-flip (our hand), `dock_volpair` to `emd_74131_m160_docked.mrc`
  (cc 0.987-0.999); `maps160/docked/conf##_emd#####_m160.mrc`, `maps160/toy_consensus.mrc` (image-weighted mean,
  corr 0.909 with conformer 1). Map scales are consistent (envelope mean density 0.068-0.074): no normalisation needed.
- Compositional variants (`make_variants.py` -> `maps160/variants/`, `rnasep_conformers10c.csv`): mobile domain =
  largest connected high-variance component (20.6 % of the molecular volume, 60 A from the centroid, 37 % of the
  variance); conf61 = conformer 49 minus that domain (8.4 % of the mass, 12 % of the conformer particles), conf62 =
  conformer 46 plus a surface blob 132 A from the domain (2.6 % mass, 7.6 % of particles).
- Real properties (`real_props.py`, |CTF| on the phase-flipped real stacks): per-particle corr with the consensus
  projection 0.404 / 0.397 / 0.157 / 0.082 / 0.051 at 117-30 / 30-20 / 20-12 / 12-8 / 8-6 A; defocus 1.29 +- 0.35 um;
  tilt 93 % within 120-150 deg.
- Simulation (`run_toy55_sim.sbatch` = toy10, `run_toy10c_sim.sbatch` = parametrised copy: SIMDIR simc, PROJ projc,
  EXTRA_MAPS variants, TAG c): N=200k, poses + CTF drawn from the refined real rows (meta_all.npz flags==1), labels ~
  image counts, 10 % blob + 10 % broken junk, `simulate_particles ctf=yes` with bfac {0,100,200} x snr {3e-4 .. 1e-2}
  chosen by `calib_snr.py` (SIGNED CTF for ctf=yes stacks; `--abs` only for flipped stacks) against the real profile;
  shuffled merge (`merge55.py`: all.mrcs, all.txt truth state=conformer, gt_particles.txt = particle conformer closed
  e1 e2); import with poses perturbed 10 deg / 3 px, refine3D refine=prob lp 8 maxits 4 from the toy consensus; the
  chain ends with `sbatch run_toy10_battery.sbatch` (env SIMDIR/PROJ/CSV/TAG propagate).

Traps hit: (1) the first sim died on `set -- $best` because calib_snr's BANDCORR line carries a "(conv .. n ..)"
suffix -> `cut -d" " -f1-6`; (2) |CTF| on a ctf=yes stack anti-correlates past the first zero (negative band 3);
(3) snr 0.02-0.1 gave low-res band corr 0.94-0.98 vs the real 0.40: SIMPLE's snr is full-image, the grid must go to
~1e-3; (4) sbatch nested in jobs inherits the environment, so plain env vars (no commas) carry parameters; (5) simulate_particles
needs an EVEN nptcls (build_refspiral throws on odd) and the oritab must end in .txt (the chain drops one particle of odd
classes into oris_in_<name>_even.txt); (6) merge55.py used to skip missing stacks silently: toy10c ran once as an
incomplete 9-conformer toy (projc_incomplete_9conf) before the merge was made to fail hard; stale all.mrcs must be
removed before a rerun because the merge is guarded by its existence; (7) refined poses after 4 refine3D iterations
from a 10 deg start: median 4.0 deg, p90 7.4 deg, shifts 1.4 px (toy10c incomplete run) = ~7 A at the periphery, the
same order as the conformational displacements, which is why per-window pose refinement is in the latent loop.

**How to apply:** every latent-model test on RNase P runs first on toy10 (and toy10c for composition) with truth:
per-particle ordering vs the conformer axis, conformer separability, closed-minority recall, window purity and
map FSC vs the truth conformers, per band. Related: [[routes-2d-to-3d-states-tests-2026-10-14]].
