---
name: klem-synthetic-benchmarks-arms-tethered
description: synthetic benchmarks for regimes beyond RNase P (user 2026-10-16): "arms" (IgM-like, 4 hinged arms with independent random swings) and "tethered" (PDH-like, 6 peripheral subunits at 50 % occupancy + 12 A tether jitter), built from REAL RNase P density pieces (not smooth shapes: user rejected the low-res first version); scripts, truth files, and the simulation driver
metadata:
  type: project
---

Generator toy55/make_synth_maps.py MODE OUT CSV TRUTH CONS (box 160, 1.464 A, 100 configurations, 1000 images each
in the CSV): core = the RNase P toy consensus recentred on its density centroid and cut by a 56 A sphere; arm/subunit
piece = the densest structured region 40-60 A from the centre (26x12x12 A ellipsoid for arms, 15 A sphere for
subunits), rigidly rotated/placed (affine_transform); arms hinge on the density surface along +-x/+-y (radii
30/41/35/13 A), theta in-plane +-40 deg, phi out-of-plane +-30 deg; subunits at surface+30 A on +-x/+-y/+-z, occupancy
0.5, displacement N(0, 12/sqrt3 A). Maps: routes/toy55/maps160/synth_{arms,tethered}/conf%02d_synth.mrc +
consensus_synth.mrc; CSV synth_{mode}.csv; truth truth_{mode}.txt (per conformer: arm angles / site occupancy+dxyz).
Driver run_toysynth_sim.sbatch MODE=arms|tethered (toy pipeline: real poses/CTF via make_oris55, bfac/snr from the ah
calibration, junk blob/broken, merge55, project proj$MODE, refine3D from the synthetic consensus with perturbed
poses; truth sim$MODE/all.txt, gt_particles.txt). Simulations trigger when f692d finishes (watcher), then the chain
(chain6/7 pattern with SRC/METAALL) and truth readouts per region (arm angle vs region coordinate; occupancy per site).
Preview figs/synth_maps_preview_v3.png. Related: [[rnasep-toy-ah-apo-holo-design]], [[klem-f692-fortran-chain]].

Trap (2026-10-17): `prg=abinitio2D` exists only in the main build ~/SIMPLE/build (2026-09-30) and ~/SIMPLE_flex2D/build, not in the kernel snapshot, klem2, stage1, coupled or kernel builds; cls_expansion runs from ~/SIMPLE_stage1/build via clsexp_stage1/run_kernel_cls.sh. sbatch --export splits commas (pass radii as 26x20x20).
Trap (2026-10-17 18:55): a project written by the kernel-snapshot build is unreadable by the older main build ("unsupported particle record width", simple_binoris); for the 2D route the synthetic projects are re-imported from the stacks with the main build (synth_2droute.sbatch: new_project + import_particles).

**Heavier sets, chain result (2026-10-18 00:xx; synth_arms2 / synth_tethered2, DEV preset, per-arm / per-site
readouts synth_perarm.sbatch):** tethered2 (subunit radius 24 A): the klem tree now localises sites: leaf states with
site occupancy 0.85-0.90 (present) / 0.27-0.28 (absent) vs ~0.5 overall, boards correlate with site 3 at 0.49 (light
set: nothing). arms2 (arm radii 26x20x20): weak: leaf std ratios 0.70-0.80 per arm, lattice-angle corr 0.47-0.57
(light set 0.74 / <= 0.48; 2D sub-classes on the light set 0.46). Verdict: composition at ~30 kDa is in reach of the
boards; hinged arms stay the 2D route's regime. Configuration purity (synth_truth.sbatch) is the wrong metric
(0.04-0.06 against 100 configurations) - judge per arm / per site.

**arms3 = Fab-mass arms (radii 37x28x28 A, ~100 kDa per arm, core 36 A; 2026-10-18):** chain leaves localise arms at
std ratio 0.66-0.69 (arms2 0.70-0.80, arms 0.74), lattice-angle corr 0.48: the per-particle tree stays weak on hinged
arms even at IgM Fab mass. Direction-bin sub-classes on arms3: best 0.55. Seeded (2D-contrast) region boards vs random
launched on its regions 1-3 (jobs 53984373/5/6). Paper check: the authors' per-conformer maps are each from the
conformer's own quality-selected subset (~50 best/mic of 21.5 M picks), rigid-core alignment, >= 50k particles
needed for 3 A by their EMCrafter test; no cross-conformer sharing.
