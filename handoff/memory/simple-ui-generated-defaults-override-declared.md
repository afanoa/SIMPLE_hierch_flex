---
name: simple-ui-generated-defaults-override-declared
description: the default value written in a UI add_input call is NOT what the built UI shows — scripts/default_audit.py generates simple_ui_default_values.f90 and add_input overwrites the declared default with it unless preserve_default=.true.
metadata:
  node_type: memory
  type: project
  originSessionId: 08367194-9f7f-4a28-aa30-3a6ae625cd31
  modified: 2026-10-01T15:13:12.704Z
---

Found 2026-10-01 fixing flex_pca's UI/commander default mismatch (Hans's email). At build time
`scripts/default_audit.py` writes `<build>/lib/simple/simple_ui_default_values.f90`
(`get_ui_default`). `ui_program%add_input` calls `apply_generated_default`, and
`ui_param%set_generated_default` OVERWRITES the declared `rval_default`/`cval_default` whenever
the table has the key. The only opt-out is `preserve_default=.true.` on the add_input.

The generator only resolves commander defaults of the form
`if( .not.cline%defined('k') ) call cline%set('k', v)` written directly in the execute path
(plus helpers listed in `TRACED_DEFAULT_HELPERS`). flex_pca sets its defaults inside
`apply_flex_pca_defaults`, which is NOT traced. So the built UI showed the `parameters`
baseline: n_probe_iters 0, maxits_pcg 2, rec_backend gridding, neigs 0, npreimages 0, even
though the source declared 5/20/gridding/10/4 and the commander used 4/4/pcg/10/4.

**Why:** reading the add_input default in source tells you neither what NICE displays nor what
runs. Three numbers can disagree: declared, generated and commander.

**How to apply:** for a commander-owned default, set the declared value AND
`preserve_default=.true.`, or trace the helper in the generator. To see what will be displayed,
run `python3 scripts/default_audit.py --source-root . --output-fortran /tmp/x.f90` and grep the
key. PUSHED to origin/master as ab6e8a149 'fixing flex UI defaults' (2026-10-01): pins n_probe_iters,
the PCG trio, npreimages, lp (16) and smpd_target (developer-only, no default); neigs still shows the
baseline 0. Same commit: state_placement=kcenter|equal_occ, preimage_auto ceiling 32 -> 8
(AUTO_NSTATES; the population-floor cap is a separate POP_FLOOR_MAX_NSTATES=32 so refine3D_states'
flex init is unchanged). Verified: fast gate 13/13 both builds, seeded A/B default + floor arms
byte-identical except the manifest state_placement label.

**NICE consequence (job builder template):** a number/text field shows the default only as a
placeholder and a blank one is NOT sent, so a wrong numeric default is cosmetic. A select
('multi') or toggle ('binary') ALWAYS submits its displayed value, so a wrong choice default
silently changes NICE-launched runs: flex_pca from NICE ran rec_backend=gridding while the
command line ran pcg, until this branch.

A trial of the generator following commander-local helpers recursively (plus commander
constants) resolved all five flex_pca values and added 21 display defaults to abinitio2D,
classify3D_refs, refine3D_auto and refine3D_states, replacing none. The trial ignores if-blocks
inside helpers, so each new entry needs checking. Not applied; the user was to raise it with Hans.
