# Why the hetero_tree root focus mask covers the whole sphere on IgM (the cheap-children offset), 2026-11-05

Side session `cheap-offset` of the IgM fixture work (coordinator: fixture session). Code read: `src_igm`, the build_igm4
source (`src/main/strategies/search/simple_hetero_tree.f90`, flex under `src/main/flex`). Work dir on batch:
`/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/cheap_offset/` (`scripts/`, `logs/`, `out/`, `dumps/`, `env/`).
Nothing under `fixture_igm/scripts/` was edited (`proto_regions.py` was read, not run), no Fortran was changed or
built, nothing was committed. Jobs: 54974143 (particle stats), 54974144 (16-node mask array), 54974825 (latent and
imprint), 54975564-66 (flex envelope), 54977541 (comparison); about 28 core-h in all.

## Verdict

- **Mechanism: a noise sort, not a background sort.** At an IgM root the ten cheap children are a noise sort.
  - They are not ordered by the particles' background level: the in-mask DC explains R2 .02-.05 of the latent.
  - flex solves its basis on the whole 600 A mskdiam sphere (92x the molecule on the fixture, 120x on real tb100k),
    so the latent's leading directions are solvent-noise modes. On the fixture the 10-dim latent carries no arm
    information (R2 <= .001), and the null twin reproduces every component (R2 .95-.99).
  - Sorting particles along PC1 imprints that direction's own basis volume into the children: a ball that fills the
    mask sphere (98 % of the difference's energy in solvent; corr with the basis volume .64 fixture, .80 real).
  - `half_variance` (3754-3759) fits one scale per map and no offset. The ball's offset (+/-.2 of the support mean,
    monotonic over the children, cross-half r 1.000) survives as offset^2 everywhere in the sphere, so the
    0.2-of-peak core is the sphere.
- **Map-side fix, per-map offset removal in `half_variance`: NO-GO.**
  - It is input-free and leaves every structural gate node unchanged (RNase P closed lineage, qtoy2, OLE: masks Dice
    >= .86 with unchanged coverage, sites identical).
  - But it only turns an uninformative 100 % mask into an equally uninformative core-only one (the 683k r1 look):
    fixture arms .87-1.00 -> .00-.06; real tb100k stays at 99.4 %.
  - It drops arm 2 on the synthetic arms root (1.00 -> .05) and changes the nulls' site discovery (qtoynull 20 -> 7).
  - Low-shell cuts, DC removal, an affine fit and radial-profile removal do no better.
- **Source-side fix, flex with an input-free envelope (`pcg_mskfile` = the node mean's own soft support): works on the
  symptom, no IgM gain, not ported.**
  - It removes the imprint at its source: offsets +/-.02, the stock rule's core 83-91 % inside the molecule, masks
    60 / 67 / 71 % of the molecule (fix / null / real) instead of 100 %.
  - But the latent still holds almost no arm information (R2 <= .024), so the mask is core-only (arms .00) as well.
  - It would change flex at every root and needs a full gate battery that the IgM evidence does not justify.
- **The whole-sphere mask is a symptom of a root latent that cannot see 50 kDa arms; the mask is not the blocker.**

## 1. Code path (line numbers: `src_igm/src/main/strategies/search/simple_hetero_tree.f90`)

- `advance`, case `S_FLEX` (613-635). A project root (`depth == 0 .and. .not. l_subtree`, 614) takes its mask from
  the 8 flex STATE maps (617-623). Every run given `infile=` is a subtree run (tb_tree, igm_tree, the 683k run, all
  fixture trees); its root takes `cheap_children` (626) and reconstruct3D of 10 states, halves only (629-634).
  Case `S_CHEAP` (636-648) puts their even/odd `_unfil` maps through `falloff_mask` (647).
- `cheap_children` (3582-3689); constants at 33 (NCHILD 10, PERBIN 1300, NBIN_MIN 20, NKM_ITS 40).
  - Latent: z1..z10 from `flex/flex_pca_coordinates.txt` (3594-3615; columns particle eo label residual
    mean_residual z1..zq).
  - Bins: view vectors from the node's ptcl3D eulers (3616-3622); 40 k-means iterations into
    `nb = max(20, nint(n/1300))` bins, 77 at 100k (3624-3657).
  - Children: the latent is centred per bin and globally (3658-3671), PC1 is the top eigenvector of its covariance
    (3672-3675), and within-bin deciles of PC1 give children 1..10 (3676-3685). No particle statistic other than the
    flex latent enters.
- `falloff_mask` (3691-3725): `core = vs > 0.2 max(vs)` (3704, MFRAC), dilation of 3 px (MDIL), Gaussian falloff
  of 6 px (MFALL), 0 outside the sphere (3708-3712).
- `half_variance` (3727-3800); its other users are `discover_sites` (3837), `variance_spread` (4376), and through
  `falloff_mask` also `oneflex_inputs` (3545), `fs_refine` (2034) and `fs_prerefine` (2201).
  - Each half map is low-passed by `fourier_filter` mode 1 (3789-3798; 4532-4559, a Gaussian with half-maximum at
    15 A whose factor at f = 0 is 1, so the DC passes untouched).
  - `sph` is the mskdiam sphere (3750) and `sup = 0.5(em+om) > 0.1 max` (3753, MSUP).
  - **The only normalisation is one LSQ scale per map fitted inside `sup` (3754-3759); there is no offset term.**
    Then `v = sum_j w_j (e_j - em)(o_j - om)` (3762-3765), smoothed with sigma 1.5 px (3766) and set to 0 outside
    the sphere (3767).
- Particle preparation, for flex_pca (`simple_flex_pca_plane_cache.f90:390/417`) and for reconstruct3D
  (`simple_matcher_3Drec.f90:440`, `prep_imgs4rec`): `norm_noise_taper_edge_pad_fft`
  (`src/main/image/simple_image_fft.f90:308`) tapers the edges, then subtracts the mean and divides by the sd of the
  pixels OUTSIDE the mskdiam disc (`lmsk`).
- flex_pca:
  - Fourier band shells 1..kmax with DC excluded (`simple_flex_pca_embed.f90:602-614`, `projected_model_kfromto`).
  - Basis solved on the whole mskdiam sphere when `pcg_mskfile` is unset. The root's flex log says so:
    "rec_backend=pcg without pcg_mskfile: the basis is solved on the spherical mskdiam support (pass
    pcg_mskfile=<solvent mask> for the envelope)". The tree's `flex_keys` (4494-4509) never sets `pcg_mskfile`.
  - Every M-step deflates 4 consensus shells, a flat-in-mask background template and a dilation template from the
    basis (`simple_flex_probe_fit_update.f90:206-290`; the root's log shows DEFLATE_BG, rank 6, removing 63-80 % of
    the basis energy).

## 2. Mechanism, with the evidence

Numpy ports matched the tree on all three IgM roots: corr(vs, tree `mask_variance.mrc`) 1.0000, Dice against the
tree's `mask.mrc` 1.000.

**(a) The children are NOT ordered by the particles' background level.** I took 8,000 particles per root and
normalised each as SIMPLE does (mean and sd of the pixels outside the 300 A disc). `din` is the normalised particle's
mean inside the disc, its in-mask DC.

| root | sd(din) | corr(din, PC1) | corr(central-disc mean, PC1) | R2 of din on the 10-dim latent | mean din per child, 1 -> 10 |
|---|---|---|---|---|---|
| fix | 0.022 | -0.12 | -0.21 | 0.016 | +0.013 ... +0.003 |
| null | 0.022 | -0.12 | -0.21 | 0.017 | +0.013 ... +0.004 |
| real tb100k | 0.033 | +0.22 | +0.33 | 0.053 | -0.019, -0.006 ... +0.004 |

The "background ordering" is a property of the child MAPS, not of the particles.

**(b) On the fixture the root latent is pure nuisance: it carries no arm information.** These are fix and null
twins: same rows, poses, CTF, quality and noise draws; the arms differ (discrete vs continuous).

| | fix | null |
|---|---|---|
| bin-centred PC shares | .30 .13 .11 .09 .08 ... | .31 .13 .11 .09 .08 ... |
| eta2 of every PC (1-10) by every arm state, or abs corr with every arm angle | <= .001 / <= .007 | <= .001 / <= .017 |
| R2 of each arm's state or angle on all 10 coordinates | .000-.001 | .000-.001 |
| twin reproduces PC1..PC10 (R2 on the twin's 10 coordinates) | .99 .97 .97 .97 .97 .96 .97 .97 .95 .96 | .99 .97 .97 .97 .97 .96 .97 .97 .96 .96 |
| corr(PC1, PC1 of the twin); cheap labels identical / within 1 | .993; 76 % / 99.4 % | |
| corr(PC1, flex contrast); flex contrast vs true amplitude Q | -.59; +.40 | -.59; +.40 |
| corr(PC1, Q / B / defocus / astigmatism / shift / flex residuals) | +.03 / .00 / +.01 / .00 / .00 / +.06-.08 | same |
| contrast projected out of the latent: new PC1 twin R2, arm NMI of its deciles | .675, .000 | .670, .000 |

PC1 follows the noise in flex's per-particle contrast estimate, and nothing in the latent tracks a 50 kDa arm. On
real tb100k, PC1 has no relation to any recorded covariate (R2 .007; contrast +.03); the PC spectrum is flatter
(.21 .12 .11 ...).

**(c) The children's maps differ by a ball that fills the mask sphere.** D is the slope of the 15 A child maps on
the child index. Values below are 10 x D / (support mean), i.e. the change over the children's range.

| root | r = 0 / 65 / 130 / 196 / 261 / 294 A, past the sphere | D energy: support / sphere solvent / outside | D cross-half r |
|---|---|---|---|
| fix | -.41 / -.63 / -.60 / -.53 / -.45 / -.30, 0 | .02 / .98 / .01 | .69 |
| null | -.47 / -.62 / -.60 / -.52 / -.44 / -.29, 0 | .02 / .98 / .01 | .69 |
| real | +.32 / +.57 / +.44 / +.39 / +.33 / +.21, 0 | .02 / .98 / .01 | .62 |
| arms root (synthetic, 110 A sphere) | -.67 / -.80 (58 A) / -.63 (102 A), -.08 at 115 A | .06 / .93 / .02 | .96 |
| r131 root (real RNase P) | +.38 (0) / +.04 flat beyond 30 A | .11 / .89 / .00 | .80 (small, not a ball) |

**(d) D is the imprint of the cheap direction's own flex basis volume.** W_y = sum_k u_k PC_k uses u, the PC1
loadings of the bin-centred latent, and the polished basis that the coordinates refer to (box 108 Fourier-padded to
200). corr(W_y, D) is .64 / .63 / .80 (fix / null / real) over the sphere and .85 / .83 / .90 inside the support.
W_y itself holds only 7 / 6 / 3 % of its energy inside the molecule support. The other basis sets (pc, merged, fitB)
correlate at |r| <= .44.

How the imprint forms:
- Particles are sorted by a noise-dominated latent coordinate, so each child's average image noise is biased along
  that coordinate's matched filter.
- Reconstruction turns the bias into +/- c_j W_y. It reproduces in both halves because both halves are selected by
  the same rule; selection makes it signal.
- The basis is solved on the whole 600 A sphere: 3.23 M voxels against a 35k-voxel molecule (92x; real 120x). Its
  noise-dominated directions are therefore solvent modes, as on the AFM RNA data (memory
  `flex-pca-solvent-mask-dominates-rna-basis`).
- The flat-in-mask deflation removes only the exactly flat component, not this radial ball.

**(e) Why `half_variance` turns that into a whole-sphere mask.**
- The far-solvent offsets per child (sphere voxels more than MDIL + 2 MFALL = 15 px from the support), as fractions
  of the support mean:
  - fix +.226 +.151 +.103 +.059 +.025 -.015 -.062 -.112 -.154 -.213;
  - null +.225 ... -.212;
  - real -.204 -.106 -.070 -.045 -.021 +.002 +.023 +.051 +.086 +.159.
  All are monotonic in the child index and identical in the two halves (cross-half r of the deviations 1.000).
- A scale fit cannot remove an additive pattern, so v is about sum_j w_j c_j^2 W_y^2. The radial profile explains
  93 / 94 / 80 % of the stock core's variance.
- 0.2 x the peak (r about 61 A fix, 35 A real) lies below offset^2 everywhere in the sphere. The core is 1.01 M /
  1.11 M / 0.39 M voxels, of which 0.3 / 0.3 / 1.1 % lie inside the molecule; the mask is 99.7 / 100 / 100 % of the
  molecule.
- The recB0 maps (after stA0) carry small offsets (+/-.05, cross-half .95-.99) that change nothing. The imprint is a
  property of latent-sorted children, not of recB0.

**(f) Where else.** The same signature appears, at varying strength, in every map set sorted by a nuisance latent:
- the synthetic arms root's cheap children: offsets +.330 ... -.344, cross-half 1.000, radial share .80;
- the nulls' flex STATE maps at their project roots:
  - qtoynull -.32 ... +.24, cross-half .999, radial share .72;
  - micnull -.12 ... +.18, cross-half .985, radial share .56;
  - olenull -.06 ... +.04, radial share .24.

It is absent from the roots with strong real heterogeneity, which dominates the latent's leading directions there:
- r131 seeds 1 and 2: -.06 ... -.01, not monotonic, cross-half .71, radial share .03;
- qtoy2 and OLE: flex-state offsets within +/-.08, radial share .01 / .07.

## 3. Map-side prototype: offset removal in `half_variance` (offset_masks.py, imprint.py)

Variants of the maps entering `half_variance`; nothing else changes (scale fit, 0.2 core, dilation, falloff):
- `solv`: each half map minus its mean over the far solvent (the `--bg yes` of `proto_regions.py`);
- `solv0`: minus its mean over the sphere outside `sup`;
- `dc`: minus its box mean (F(0) = 0 in the 15 A filter);
- `aff`: scale and offset fitted together over the sphere;
- `hp1/2/3`: shells with q < c / mskdiam removed;
- `rad`: minus its radial profile inside the sphere.

Root masks (mask > 0.5 as % of the molecule support; coverage = fraction of each true envelope inside the mask:
`vols/core.mrc` > 0.1 max, `arm_masks.mrc` 1-5):

| variant | fix: mask / core / a1 a2 a3 a4 a5 | null: mask / core / a1-a5 | real: mask | core voxels inside the support (fix / null / real) |
|---|---|---|---|---|
| stock | 99.7 / 1.00 / .87 1.00 .99 .91 .97 | 100 / 1.00 / .92 .99 .99 .92 .96 | 100 | 0.3 / 0.3 / 1.1 % |
| solv | 85.4 / .96 / .00 .04 .00 .06 .00 | 90.7 / 1.00 / .00 .07 .00 .20 .09 | 99.4 | 13 / 4 / 13 % |
| solv0 | 83.6 / .95 / .00 .03 .00 .00 .00 | 90.1 / .99 / .00 .05 .00 .20 .09 | 99.4 | |
| dc | 96.9 / 1.00 / .45 .43 .42 .43 .50 | 99.2 / 1.00 / .57 .57 .45 .46 .55 | 100 | |
| aff | 57.6 / .73 / .00 x5 | 57.4 / .73 / .00 x5 | 98.5 | |
| hp1 | 77.1 / .86 / .07 .87 .18 .33 .17 | 86.5 / .90 / .50 .98 .73 .82 .76 | 93.6 | 9 / 2 / 35 % |
| hp2 | 80.3 / .91 / .00 .03 .00 .00 .02 | 83.0 / .93 / .02 .06 .00 .03 .27 | 69.1 | 29 / 13 / 67 % |
| hp3 | 79.2 / .91 / .00 x5 | 84.2 / .95 / .01 .00 .01 .05 .06 | 88.8 | 65 / 19 / 64 % |
| rad | 59.7 / .74 / .00 x5 | 69.7 / .82 / .00 x5 | 80.6 | 88 / 78 / 47 % |

- Removing the offset does stop the whole-sphere core.
- What remains is the rest of the same imprint: a blob centred near the core, still mostly solvent. Every variant
  gives the core-only mask (the 683k r1 look) or a partial blob. None covers the arms, because the latent holds no arm
  signal (2b).
- On real tb100k the molecule coverage barely moves (99.4 % with `solv`): the imprint around the molecule still
  covers it.
- The recB0 masks and sites of the three roots are unchanged by `solv` / `solv0` / `dc` (mask Dice .97-.99, sites
  .98-1.00). `aff` moves sites (Dice .13-.73 on some), another reason to drop it.

## 4. Gates: the same before / after on saved maps (no reruns needed; every tree below kept its maps)

Nodes: `hetero_tree_20261028/math_reworks/gates/base_*`, `route9/guardfix/runs/r131_md1_s1`,
`subcls_mask/arms/tree/a`. Ports vs the trees' own masks: Dice .966-1.000, variance corr 1.0000. `solv` vs stock:

| node (mask maps) | mask, % of the molecule (Dice) | truth coverage | recB0 sites (best Dice to stock) |
|---|---|---|---|
| r131 seed 2 root (cheap) | 94.0 -> 93.9 (.876) | closed .97 -> .99, open .98 -> .98 | 4 -> 4 (.99-1.00); closed-arm site .83 -> .83 |
| r131 seed 1 root (cheap) | 93.5 -> 92.2 (.860) | closed .98 -> .98, open .97 -> .96 | (labels not kept) |
| r131 r1 (parent recB0) | 61.2 -> 61.2 (.995) | closed .98 -> .98 | 5 -> 5 (.98-.99); closed .96 -> .97 |
| r131 r2 (parent recB0) | 60.5 -> 60.5 (.996) | closed .95 -> .95 | 3 -> 3 (1.00); closed .86 -> .86 |
| qtoy2 root (flex) | 51.7 -> 51.7 (.993) | - | 6 -> 6 (.99-1.00) |
| qtoy2 r1 / r2 (parent recB0) | 90.0 -> 90.0 (.998) / 29.4 -> 29.6 (.995) | - | 14 -> 14 (.95-.99) / 4 -> 4 (.98-1.00) |
| OLE root (flex) / r1 | 45.1 -> 46.4 (.979) / 33.3 -> 33.4 (.996) | - | 4 -> 4 (identical) |
| qtoynull root (flex) | **98.6 -> 79.6 (.749)** | - | **20 -> 7 (.13-.88)** |
| olenull root (flex) | 99.7 -> 99.8 (.937) | - | **7 -> 5 (.91-.97)** |
| micnull root (flex) | 100 -> 100 (.982) | - | 0 -> 0 |
| arms root, synthetic (cheap) | **99.8 -> 45.2 (.265)** | **arms 1.00 1.00 1.00 1.00 -> .52 .05 .83 1.00** | 4 -> 4 (>= .99); arm sites .95 .94 .00 1.00 unchanged |

- RNase P closed-state lineage: unchanged (closed arm in every mask and site as before).
- qtoy2 and OLE: unchanged.
- Nulls: the mask and the site candidates change wherever their maps carry the reproducible offset (qtoynull most).
  The split decision would need replays.
- Arms root: the structural peak (arm 4) survives, but arm 2 drops below 0.2 of it and leaves the mask. With the
  stock 99.8 % mask the tree's sites find arms 1, 2 and 4.
- For `island_gates.md` ("bg yes" variants): with offsets removed, the stock `discover_sites` port on the nulls'
  recB0 maps changes too, even where no island appears.
  - qtoynull: 20 -> 7 sites (best Dice to the stock sites .13-.88).
  - olenull: 7 -> 5.
  - Any port that removes offsets before discovery changes the nulls' seed candidates, so the null replays must
    include it.

## 5. Particle-side idea (offset / DC removal of the particles before the latent): not supported

The particles' in-mask DC explains R2 .016 / .017 / .053 of the latent (2a), and flex already excludes the DC shell.
Projecting the particles' contrast out of the latent leaves a nuisance PC1 that the twin still reproduces (R2 .67)
and that carries no arm information (2b). There is no particle background to remove that would change the children.

## 6. Source-side prototype: flex with an input-free envelope (env_flex.sbatch; no Fortran change)

Setup: existing keys only, no Fortran change.
- Envelope (`env_mask.py`): the node mean (`mean/recvol_state01.mrc`) low-passed at 15 A, thresholded at MSUP x its
  maximum inside the sphere (`half_variance`'s `sup`), then made soft as `falloff_mask` makes a core (MDIL px
  dilation, MFALL px Gaussian falloff, 0 outside the sphere).
  - Support 35.4k / 35.1k / 27.0k voxels (fix / null / real) in 11 / 8 / 1 components.
  - Mask > 0.5 covers 7.9 / 7.2 / 4.9 % of the sphere; flex reports a hard support of 7.3 / 6.5 / 4.2 % of its
    108-px box.
- The root's flex_pca was rerun with the tree's own keys (from `node_r/flex/hetero_step.sh`) plus
  `pcg_mskfile=env.mrc`: 12-13 min, as stock.
- The cheap children were rebuilt by the tree's rule (`cheap_labels.py`) and reconstructed with the tree's cheap keys
  (`selection` + `reconstruct3D`), then given the stock mask rule (`offset_masks.py`).

| root | child offsets (far solvent / support mean) | cross-half r | stock-rule core: voxels (share inside the molecule) | mask, % of the molecule | coverage core / a1-a5 |
|---|---|---|---|---|---|
| fix, stock flex | +.226 ... -.213, monotonic | 1.000 | 1,007,075 (0.3 %) | 99.7 | 1.00 / .87 1.00 .99 .91 .97 |
| fix, envelope flex | -.012 ... +.014, not monotonic | .90 | 5,721 (90 %) | 60.2 | .75 / .00 x5 |
| null, stock flex | +.225 ... -.212 | 1.000 | 1,107,287 (0.3 %) | 100 | 1.00 / .92 .99 .99 .92 .96 |
| null, envelope flex | -.013 ... +.018 | .93 | 5,872 (91 %) | 67.4 | .82 / .00 x5 |
| real, stock flex | -.204 ... +.159 | 1.000 | 390,693 (1.1 %) | 100 | - |
| real, envelope flex | -.034 ... +.026 | .99 | 3,036 (83 %) | 70.8 | - |

- Every offset variant of section 3 is now inert (mask Dice >= .96 against stock): there is no offset left.
- The children still differ by the imprint of their own latent direction: D matches the cheap direction's polished
  basis volume at r .86 / .86 / .77 over the sphere and .98 / .98 / .97 inside the support. But that imprint now sits
  inside the molecule (61 % of D's energy in the support, stock 2 %): a compact mode at the box centre, D at r = 0
  +1.06 / -1.10 / +1.07 and |D| <= .07 beyond 30 A on the fixture.
- Latent with the envelope:
  - fix: PC shares .33 .20 .15 .10 .05. PC1 no longer follows flex's contrast (r -.007, stock -.59); all covariates
    give R2 .006.
  - Arms (fix): every PC's eta2 <= .009. R2 of each arm's state on all 10 coordinates .024 / .011 / .005 / .003 /
    .002 (stock <= .001). Continuous angles .000. The cheap deciles' NMI with every arm is .000.
  - null: arm R2 <= .002.
  - real: flat spectrum .16 .11 .11 .11 .10; PC1 vs covariates R2 .08 (shift -.13, contrast -.12).
  - Twin (envelope fix vs envelope null): PC1..PC10 R2 .94 .89 .82 .86 .67 .76 .63 .46 .55 .71 (stock flex .99-.95). The leading directions are still nuisance shared by the twins; only the trailing, less shared components hold the little arm signal there is (a1 eta2 .009 on PC8, twin R2 .46).
  - Basis energy inside the molecule support per component, stock vs envelope: shares of each polished component's energy inside the support (the rest is sphere solvent, which with the envelope includes its ~10 px falloff):
    - stock fix .44 .09 .16 .09 .05 .08 .07 .06 .03 .04; real .07 .04 .04 .01 .04 .03 .03 .04 .03 .02;
    - envelope fix .72 .37 .47 .22 .32 .20 .24 .31 .21 .15; real .52 .16 .14 .13 .12 .14 .18 .19 .16 .15.

Reading:
- The envelope removes the whole-sphere mask at its source: no solvent imprint, molecule-internal cores, and PC1
  stops tracking contrast noise.
- At 50 kDa per arm the root latent still holds almost no arm information. Its leading direction becomes a compact
  noise mode at the centre, and the mask ends core-only, exactly where the map-side offset removal ends.

## 7. Recommendation

1. **NO-GO: per-map offset removal in `half_variance`, or its DC / affine / low-shell / radial variants.**
   - It turns the IgM root mask from 100 % into a core-only mask, arms .00-.20 (real tb100k: 99.4 %).
   - It drops arm 2 on the arms root and changes the nulls' site discovery (qtoynull 20 -> 7, olenull 7 -> 5).
   - It leaves RNase P, qtoy2 and OLE untouched.
   - It removes a symptom, not the cause, and buys nothing for IgM.
2. **Hold, not ported: the flex envelope** (`pcg_mskfile` = the node mean's own soft support).
   - Pros: it is the clean, input-free fix of THIS symptom, with no new constants and an existing key; it works on
     fixture, null and real.
   - Cons: the resulting mask is core-only as well, so it cannot help IgM states.
   - It changes flex at every node that runs flex: project roots' state masks, fflex (whose `pcg_mskfile` is the
     node mask, `submit_fflex` 1092), discovery and seeds. A full gate battery would be needed: r131 closed lineage
     md1 both seeds, qtoy2 md1, qtoynull / olenull / micnull md0, OLE md1; a build_igm6 plus about 8 trees,
     roughly 300 core-h. The IgM evidence does not justify that under the design freeze.
   - Keep it as a candidate only if a non-IgM reason appears. The nulls' noise-sorted flex STATE maps (qtoynull
     offsets .999 reproducible, 20 site candidates) are the obvious one.
   - Risk to check if ever tested: the envelope is the 0.1 x max support plus about 10 px, so low-density mobile
     parts outside it leave the basis (memory `focus-mask-truncates-mobile-tip-falloff-prior-2026-10-21`). Real
     IgM's envelope is a single 27k-voxel component, with smeared arms.
   - The minimal port, for reference: in `advance` case `S_MEAN` (607-612), write the envelope from
     `mean/recvol_state01.mrc` (the `sup` rule plus `falloff_mask`'s soft edge), then
     `cl%set('pcg_mskfile', ...)` before `submit(i, S_FLEX, cl)`. Build in a new build_igm6 by the
     `tree_build4.sbatch` patch-and-build pattern.
3. **What the IgM root mask is telling us.** The whole-sphere mask is the visible symptom of a root latent with no arm
   signal (fixture arm R2 <= .001 stock, <= .024 with the envelope). The mask is not the blocker; arm information has
   to come from somewhere other than the root's flex latent. The fixture session's per-site focused roots found a1 at
   NMI .128, and the cross-validated stage B is the open prerequisite there.
4. Optional, log-only (no decisions, user's call): `falloff_mask` could log its maps' far-solvent offset range and
   their cross-half r, the way `variance_spread` logs at the split.
   - Large, reproducible offsets flag a noise-sorted node: IgM +/-.2 (r 1.000), the arms root +/-.33,
     qtoynull's states -.32 ... +.24.
   - Structural nodes stay small: r131 within +/-.06, qtoy2 / OLE within .08.

## Files

- Scripts (batch, `cheap_offset/scripts/`):
  - `bg_latent.py`: particle in-mask DC vs PC1, twin check.
  - `latent_cov.py`: PCs vs covariates and truth, twin R2, contrast projected out.
  - `offset_masks.py`: half_variance / falloff_mask / discover_sites port with the offset variants, any node.
  - `imprint.py`: D, W_y, filter variants.
  - `env_mask.py` and `cheap_labels.py`: envelope and cheap-children mirrors.
  - `basis_support.py`: basis energy inside the support.
  - sbatch wrappers: `bgl`, `om` (+ `om_spec.txt`), `imp`, `env_flex`, `envcmp`.
- Reports:
  - `out/bgl/*_bg_latent.txt`, `out/om/<node>/report.txt` (+ masks / variance maps for stock and `solv`);
  - `out/imp/*_imprint.txt`, `*_latent_cov.txt`;
  - `env/<ds>/` (flex run, cheap maps, `om/report.txt`, imprint, latent).
