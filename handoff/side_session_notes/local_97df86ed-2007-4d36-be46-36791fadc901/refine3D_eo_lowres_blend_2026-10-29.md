# refine3D low-resolution even/odd reference blend (tree-micelle-synth, 2026-10-29)

## What it does
- `blend_lowres_eo_for_registration` (contained in `read_mask_filter_reproject_refvols`, `src/main/strategies/search/simple_matcher_refvol_utils.f90`) replaces both half references below a cutoff with their average.
- The same blended pair feeds the polar tables and the Cartesian (refine=cont) references.
- **Cutoff:** `find4eoavg = max(K4EOAVGLB = 4, Fourier index at FREQ4EOAVG3D = 20 A)`, lowered to the shell where the state's FSC falls to `FSC4EOAVG3D = 0.95` (`src/defs/simple_defs.f90`).
  - At box 88 / 3.545 A that is shell 16 or lower, i.e. everything coarser than about 20 A.
- **Effect when it acts:** an even particle's reference holds odd-half data, and its own half's, below the cutoff. This is the usual gold-standard low-resolution averaging for registration.

## When it acts
- **Only when `lp=` is NOT on the command line.** The routine returns at once when `params%l_lpset`, and `l_lpset = cline%defined('lp')` (`simple_parameters_phases.f90:591`). solve3D stages set `l_lpset` themselves.
- It therefore does **not** act in:
  - any hetero_tree node refinement (stA0 / stage A / stage B pass lp 15 / 10);
  - the fixed-pose state classifier runs of the micelle fixture (lp=10);
  - the soft multi-state prototype (lp=10).
- **No result of the micelle fixture work depends on the blend** (but see below: with lp set there are no separate halves at all). My earlier message to the coordinator saying the soft runs were "not cross-half below 20 A" was wrong and is corrected here.

## When it matters
- **Any cross-half scheme run without an explicit lp:**
  - scoring each half against the other half's references (cross-half classification, a cross-half posterior);
  - an e/o-independence argument that needs the references to hold no own-half data below 20 A.
- **What to do:** pass `lp=`, or skip the blend. The `SIMPLE_SOFT_SHRINK` reference path in the soft2 scratch build always skips it.

## The bigger point: with lp= the matching references are the MERGED map in both halves
- **Where:** `read_mask_filter_refvols` (same file, lines ~433-446). The branch is `else if( params%l_lpset .and. .not. l_nonuniform_mode .and. .not. l_diag_halves )`: it reads `vols(s)`, "the average volume that will occupy both even and odd", and copies it to `vol_odd`. (Nonuniform mode has its own merged-reference rule, `l_use_merged_nu_ref`.)
- **Consequence:** whenever `lp=` is given, every particle is scored against a map that contains its own image, at every shell. There are no independent half references.
- **Who passes lp:**
  - every hetero_tree node refinement (stA0 / stage A / stage B, lp 15 / 10);
  - the micelle fixture's fixed-pose classifier;
  - all `SIMPLE_SOFT_STATES` runs on the soft and soft2 builds.
- **How it surfaced:** on the soft2 build, `SIMPLE_SOFT_SHRINK` logged an e/o difference FSC of exactly 1.000 at every shell and iteration. Its runs reproduced the unshrunk runs exactly.
- **soft3 scratch build:** under `SIMPLE_SOFT_STATES` with fixposes=yes, the true even/odd halves are kept even when lp is set (`l_soft_halves`). Nothing else changes.
