# Port proposal: an add-only per-site split for hetero_tree (2026-11-06)

For approval. Nothing is built, launched or committed.
- Evidence: tree_igm_fixture_report.md (BATTERY, SCORECARD UPDATE, GATES); discrete_vs_continuous.md sections 3-7b; island_gates.md.
- Python reference code (batch): fixture_igm/scripts/xsbg.py and proto_regions.py; fixture_igm/discrete/scripts/dc_battery.py.

## In one paragraph
The tree's own tests sometimes stop a node above the depth cap. The layer looks only at such nodes, and asks of each site whether the particles fall into two discrete states there. The sites are the node's own plus islands: variance peaks away from the core. For each site:
- A cross-validated two-state classification splits the particles into two groups.
- GAPw calls DISCRETE only when the consensus density between the two lobes of the group difference is nearly empty.

On a DISCRETE call, the leaf becomes a split node, and the two groups become children that run the normal node pipeline. Otherwise the leaf is filed as before. The tree's own splits, sites and seeds never change, and the layer is off by default.

## Price (the lean scope)
- **Where it runs:** only at would-be leaves above the depth cap that hold at least 2 x minpop particles.
  - **692k RNase P recipe** (f692_lf_s1: maxdepth 4, minpop 3000, 5 h 14 min wall):
    - It reaches 12 of the 66 leaves: the depth-3 leaves r412-r542, 8.7-16.6k particles each.
    - That is about 50-60 site jobs, at the usual 4-5 sites per node.
    - Every closed leaf sits at depth 4, out of the layer's reach.
  - **IgM fixture:** the stock tree's root is such a leaf, because the locality guard merges its split. So the layer runs there on 100k particles with 6 sites (the core site + 5 arm islands).
- **Cost per site:**
  - Python, measured: 5-10 min on 16 cores at 12-53k particles (gates); ~35 min and ~10 core-h at 100k (battery).
  - Fortran: expected cheaper, because the poses are fixed. The consensus projections and the noise model are computed once. Each of the 10 iterations is 4 half reconstructions at a 12 A crop, plus 2 projections per particle. Gate 1 measures it.
- **On f692:**
  - Core-h: about 80-100 core-h at the Python rate.
  - Wall: the layer runs as each eligible leaf finishes. It adds about one site job (5-10 min) to the wall time only when such a leaf is on the critical path.
- **The expensive case:** a DISCRETE call creates two children that run the full node pipeline (~1 h each). The gates found no such call on RNase P or the nulls.

## Why each part (knockouts)
All on the IgM fixture (5 discrete arms, 100k) and its null twin (the same arms, continuous).

| part | without it | report section |
|---|---|---|
| islands | The tree's sites sit only in the core: fixture root 1 site, null 5 sites, all core. No arm is ever examined. | Task 2; findings 17:03-18:20 |
| site-focused classification | Given the island sites, the existing site layer finds no pole on any of the five arms (endpoint z 0.4-1.8). | (b) + (d) replay |
| cross-validation (per-half seeds, swap, half alignment) | With merged references, the null's per-arm splits pass the split test and the guard while tracking nothing (groups 65 / 59 %). This is the noise-lock. | NULL per-arm controls |
| GAPw | Every map-level gate also passes the null's real continuum cuts. 15 A relative passes null a1 (2.2-2.6) and a5 (2.5); cross-half z passes everything; rigid-core controls split at z 26. The null ends in 4 groups, the largest ~25 %. | BATTERY |
| add-only, at leaves | As the only split rule, GAPw would stop RNase P at r21, because the closed lineage's own splits read not discrete. | GATES, wave 2 |

## What it buys
- **Fixture:** the a1 arm splits into groups 72 / 28 % pure (NMI .146; seed 2 .138; GAPw .08-.10 in both seeds). The other four arms stay unsplit, and the null twin stays one leaf.
- **RNase P and the null fixtures:** nothing. 89 gate runs gave 0 DISCRETE calls; no closed state was gained or lost.
- **Misses, by design or lack of power:**
  - three states in a row (a2, a5), because the middle state fills the gap;
  - the +-20 deg arms a3 / a4, which the classification does not find;
  - nodes under ~50k particles for a 50 kDa arm (a1 at 25k reads .32).
- **Untested:** other continuum shapes (ctun, ctpk), the three-domain mm3q, and real IgM.

## Design

### Where it acts
- In `decide`, take a node about to be filed as a single leaf. That happens when the split test fails, a group falls under minpop, or the guard or floor merges the split. If the node is at depth < maxdepth and holds >= 2 x minpop rows, its leaf filing waits for the layer.
- The layer never acts at a node the tree splits, or at the depth cap.

### Sites
- The node's own sites (sites/site_XX.mrc), as discovered.
- Islands, written to separate files (sites/island_XX.mrc, islands.txt) that only the layer reads. The tree's site ranking and seeding therefore stay byte-identical.

The island rule (proto_regions.py; the ADDITIVE reading in island_gates.md):
1. Compute the half-variance from half maps with each map's solvent mean removed. The solvent is everything inside mskdiam that lies beyond MDIL + 2 MFALL of the support.
2. Claim the stock core (above MFRAC of the maximum) first.
3. Take every variance peak inside the support, in decreasing order. Its connected piece above MFRAC of the peak becomes an island if it touches nothing claimed and holds >= DMINVOX voxels.
4. Keep an island when discover_sites' SVD factor test keeps >= 1 factor (z >= DFACT_Z, share >= DFACT_SHARE). Its mask is dilated by MDIL with an MFALL falloff, and DMAXMOL applies.

The rule adds no new constants. The solvent-mean removal is what finds the fixture's a1 island. The nulls have no islands with or without it.

### Two-state classification per site (xsbg.py, swap design)
- **Inputs:** the node's stage-B project (poses, CTF, even / odd), the site mask, and the consensus (the node's rows reconstructed at their stage-B poses, halves kept).
- **Scoring, at fixed poses and 12 A:**
  - each particle's residual to its scaled consensus projection;
  - each state's site-masked deviation M (R_k - mean), projected with CTF and shift;
  - a whitened Gaussian log-likelihood, with per-shell noise taken from the residuals;
  - hard assignment.
- **Cross-validation:**
  - each half gets its own random seed labels;
  - each half is scored against the other half's state maps;
  - before every iteration, half 1's two labels are kept or swapped, whichever better matches half 0's site deviation.
- K = 2, 10 iterations. Output: two groups and their per-half maps.

### The DISCRETE call (GAPw, discrete_vs_continuous.md section 3)
Computed inside the site, on the two groups' half maps and the consensus, Fourier-cropped to ~5.1 A/px with a 15 A low-pass:
- **G1 (a signal-dominated factor):** after nuisance terms are removed, the first cross-half factor has z1 >= 3 and r1 >= .5.
- **Lobes:** w = sign(F1e + F1o) sqrt(max(F1e F1o, 0)); M+ = w > .5 max w, M- = w < .5 min w.
- **G2 (density moves):** the lobe centroids are >= 10 A apart.
- **G3 (lobes reproducible):** the even and odd lobe axes agree within 30 deg.
- **B:** the lobes' convex hull minus the lobes dilated by 10 A. It needs >= 10 voxels, else no call.
- **GAPw:** the consensus density in B over the weaker lobe's density, per half, averaged.
- **Call:** DISCRETE when G1-G3 pass and GAPw < .28.

Two notes on the call:
- **Readout:** the port uses the two groups directly (the "direct" readout). The validated six-bin readout gave the same call on all 104 runs (14 battery, 57 + 33 gate rows). The direct one needs no extra ordering pass or reconstruction.
- **Constants:** all fixed, none a knob, all in A: tau .5, 10 A dilation, z1 >= 3, r1 >= .5, 10 A, 30 deg, 10 voxels, .28. The .28 is the fixture null's median minus 3 robust sd.

### The split
- If several sites are DISCRETE, the lowest GAPw wins: one split per step.
- The leaf becomes a node with two children (`add_child`). They run the full pipeline, and the layer again if they too stop above the cap.
- The split is logged as "site k discrete (GAPw .., z1 .., L .. A)". Per-site statistics go to sites/xs_<k>/stats.txt.

### Switch
A new hetero_tree key, `site_split = yes | no`, default no:
- With no, nothing changes; gate 0 checks byte-identical leaves and labels.
- RNase P recipes keep it off.

## Code

**Base:** fork/hetero-tree 72e5aaba0 (the f692 combined code), plus the IgM build's two functional patches (tree_stage_box, tree_mem_rate). The two diagnostic mask patches stay out.
- The work goes on a new local branch in its own worktree, and builds in a new dir on batch (build_igm5). build_igm2-4 stay untouched.
- Nothing is committed until you review it.

| file (src/main/...) | change | lines |
|---|---|---|
| strategies/search/simple_hetero_tree.f90 | Islands: a pass after discover_sites' region loop (:3902), written apart. The key. A node step S_SITESPLIT (21): `decide` (:2786) defers the leaf and submits one worker per site, then adds two children or files the leaf. | ~250 |
| strategies/search/simple_site_split.f90 (new) | `run_site_split(cline)`: reads the node project and the site mask, runs the classification and GAPw (with a small 3D convex hull for the voxel test), and writes labels and stats only. It follows the layering policy: the worker never writes the project. | ~800 |
| commanders/simple/simple_commanders_refine3D.f90 | `commander_site_split`, a thin execute(cline), like commander_site_discover (:77, :3060) | ~25 |
| exec/simple_exec_refine3D.f90 | case('site_split'), next to site_discover (:97) | ~3 |
| ui/simple/simple_ui_heterogeneity.f90; params | the key, and site_split as an internal program like site_discover (:852) | ~35 |
| tests | unit tests: label alignment, the convex-hull voxel test, GAPw on synthetic two-lobe vs continuum maps | ~200 |
| doc/implementation_notes/hetero_tree_site_split_2026_11.md | the method and the gates | - |

That is about 1,300 lines, two thirds of them in the new worker. It is a real addition to a tree you already find complicated, but it is confined to leaves and off by default.

## Acceptance (Fortran, build_igm5)
Staged: gates 0-1 run after approval; gates 2-6 run only after you have seen 0-1.

| # | run | passes when | ~core-h |
|---|---|---|---|
| 0 | site_split=no vs the base build, on qtoy2 (md1) | byte-identical leaves and labels | 25 |
| 1 | the worker alone, on the saved gate preps (89 runs) and on the battery's fixture / null a1 and a5 | Python's calls: 0 DISCRETE on the gates, fixture a1 DISCRETE in both seeds, null no call. Group agreement with the Python labels >= .9 at the strong factors (r211 s01, r21 s01-03). GAPw within .05. Also measures the cost per site. | 20 |
| 2 | IgM fixture and null twin, full trees, site_split=yes | the fixture splits a1 and no other arm; the null stays one leaf | 110 |
| 3 | qtoynull, olenull, micnull trees | one leaf each (largest >= 80 %) | ~100 |
| 4 | r131 corner trees (bc_r131, both seeds) | corner as in the base build | ~50 |
| 5 | first time anywhere: ctun, ctpk, mm3q (60k roots) | ctun: no DISCRETE; mm3q: DISCRETE at its motion sites; ctpk: reported | ~130 |
| 6 | optional: the f692 recipe with site_split=yes | closed leaves identical; any added split only at the 12 depth-3 leaves; wall <= 15 min over 5 h 14 min | production-sized |

Gates 0-5 total ~435 core-h. A failed gate stops the sequence, and no constant is tuned to pass one.

## Risks
- **Value:** one arm of five on one fixture, nothing on RNase P. The port pays off only if mm3q or real IgM-like data give calls.
- **GAPw's threshold** comes from one null shape, a fan continuum. ctun and ctpk (gate 5) are its first test elsewhere.
- **Real IgM:** about 1/3 of the 683k consensus particles sit in the wrong 60-degree register. That nuisance is itself a two-state contrast at the J slot and the arm sites, so any DISCRETE call there must be checked against register before it counts.
- **Independent sites:** discrete sites split one per level, so each later site sees half the particles. A product over sites at one leaf (mm3q's lesson) is not proposed: no data have yet shown two calls at one leaf.
- **Island variant:** the nulls were gated with the solvent-mean removal (no islands either way). The RNase P gate used islands without it. On f692 the closed leaves are out of reach regardless.

## Not proposed
- **Running the layer at every node** (a site product over the tree): more power, but it would act where the tree already splits. The 692k root and r5 were never gated.
- **A Python post-pass on finished trees with kept maps**, as the gates ran: no Fortran, but it would report DISCRETE leaves for a subtree rerun rather than split them in place. This is the fallback if the port's size is not justified.

## Decision
Approve the design and the build with gates 0-1. The design: leaf-only, islands for the layer only, the direct GAPw readout, off by default.
