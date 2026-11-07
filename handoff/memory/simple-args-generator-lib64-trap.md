---
name: simple-args-generator-lib64-trap
description: "New parameters keys are silently rejected by workers ('X argument is not allowed') when the build dir has lib64/: simple_args_generator.pl writes there while CMake compiles lib/simple/simple_args.f90; copy or regenerate into lib/ before make"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-07T21:14:44.225Z
---

Found 2026-09-07 adding `nfits` and `pindfile` to `simple_parameters.f90` for the flex_pca
worker job_descr (Step C of the distributed alignment plan). The command-line checker
(`simple_cmdline :: parse_command_line_value`, `allowed_args%is_present`) uses the GENERATED
`simple_args.f90`, produced by `scripts/simple_args_generator.pl` from the parameters file.
CMake's custom command declares the output as `<build>/lib/simple/simple_args.f90`, but the
perl script writes to `<build>/lib64/simple/` whenever `lib64` exists (it does after the first
install on the box), so from then on the compiled list is frozen at its first generation and
every later key is "not allowed" for distributed workers while shared-memory runs are fine
(they never parse the key). Registration itself needs only the params field plus
`reg%add_int/add_file` in `simple_parameters_parse.f90`; no UI entry (pcafit has none).

**How to apply:** after editing `simple_parameters.f90`, run the generator and copy
`lib64/simple/simple_args.f90` over `lib/simple/simple_args.f90` before `make` (the box
`~/gate_step.sh` does this); verify with `grep <key> build/lib/simple/simple_args.f90`.
Also: new source files need `cmake .` before make (no CONFIGURE_DEPENDS on the source globs).
Related: [[flex-pca-em-split-and-syntax-check-recipe]].
