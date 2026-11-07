# hetero_tree: the add-only per-site split (`site_split=yes`), 2026-11-06

Ported from the IgM fixture prototypes (`xsbg.py`, `proto_regions.py`, `dc_battery.py`; the proposal and its evidence in the
session notes `port_proposal_site_split.md`, `tree_igm_fixture_report.md`, `discrete_vs_continuous.md`, `island_gates.md`).
Off by default; with `site_split=no` the tree is unchanged.

## Where it acts

`decide` (`simple_hetero_tree.f90`): a node the tree would file as one leaf (the split test merged everything, `ng == 1`), at a
depth below `maxdepth` and holding at least `2 x minpop` non-junk particles, is not filed. One `site_split` worker runs on its
stage-B project (step `S_SITESPLIT`, directory `xsplit/`). When it returns, `decide` reads `xsplit/summary.txt`: a DISCRETE site
(the lowest GAPw of several) whose two groups both hold `minpop` particles splits the node into two children that run the full
pipeline (and the layer again if they stop above the cap); otherwise the leaf is filed as before. Junk states become junk leaves
either way. Nodes the tree splits, and nodes at the cap, never see the layer.

## Sites

The node's own sites (`sites/site_NN.mrc`) and the islands (`sites/island_NN.mrc`, `islands.txt`), which only this layer reads:
`discover_islands` recomputes the half-map variance with every map's solvent mean removed (`half_variance(l_bg=.true.)`: the solvent
is the sphere beyond `MDIL + 2 MFALL` voxels from the support), claims the stock core, then takes every variance peak inside the
support in decreasing order; its connected piece above `MFRAC` of the peak is a region when it touches nothing claimed and holds
`DMINVOX` voxels, and an island site when `discover_sites`' SVD factor test keeps a factor (`DFACT_Z`, `DFACT_SHARE`) and the
mask (dilated by `MDIL`, `MFALL` falloff) covers at most `DMAXMOL` of the molecule. No new constant.

## The worker (`simple_site_split.f90`, program `site_split`)

Per site, on `state_mergeup`: `site_focus` (the node mask times the site mask), then `site_split_em`: per-half random starts, ten
iterations of the cross-half EM on the site path (`cross_em(l_site=.true.)`, the full band): each half's references are its
consensus plus the site-masked deviation of each group from the groups' weighted mean (`site_refs`), so a particle is scored on
its consensus residual against the other half's site deviations, as `xsbg.py` does; the prior is flat and the labels hard (soft
weights at a noise-level site contrast pull every group map onto the consensus, and the standard mixing prior then drains the
smaller group within three iterations); half B's groups are renumbered to half A's by the cross-half correlation of their
deviations before every E-step; the group maps are neither docked onto the other half's consensus nor FSC-shrunk (a shift per group
map makes the two references a consensus-gradient dipole, and the EM then converges to a translation split that both halves
reproduce). The two groups' unmasked half maps come from the final labels.

The call (`gap_call`), on the 15 A low-passed group half maps inside the site:
- G1: the first cross-half factor of the site deviations is signal-dominated after the nuisance terms (the mean, the consensus,
  its Laplacian and rotational blur projected out): z1 >= 3 and r1 = tanh(z1 / sqrt(n_eff - 3)) >= 0.5.
- Lobes: w = sign(F1e + F1o) sqrt(max(F1e F1o, 0)) with F1 the group difference; M+ = w > 0.5 max w, M- = w < 0.5 min w.
- G2: the lobe centroids are at least 10 A apart (half-separation); G3: the even and odd lobe axes agree within 30 deg.
- B: the convex hull of the lobes (`hull_mask`, an incremental 3-D hull of the lobes' boundary voxels) minus the lobes dilated by
  10 A; at least 10 voxels. GAPw = the consensus density in B over the lower lobe density, per half, averaged.
- DISCRETE when G1-G3 and B pass and GAPw < 0.28. Constants from the prototype, none a knob.

Outputs in `xsplit/`: `summary.txt` (`site n1 n2 z1 r1 lsep_A axis_deg nB GAPw GAPw_even GAPw_odd call fails`) and per site
`<tag>/labels.txt` (every project row, 0 = none) and `stats.txt`. The project is never written.

## Acceptance

The proposal's gates: 0 (`site_split=no` byte-identical leaves and labels on the gate fixtures), 1 (the worker alone against the
Python calls on the saved gate preps and the battery), then 2-6 (IgM fixture and null twin, the nulls, r131, ctun/ctpk/mm3q,
optionally f692). Results are recorded here as they come in.
