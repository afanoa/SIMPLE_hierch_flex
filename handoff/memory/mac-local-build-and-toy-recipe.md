---
name: mac-local-build-and-toy-recipe
description: when the cluster is unreachable the user says "test locally" - Mac build of SIMPLE_clsrender needs the macOS 26.5 SDK (27.0 breaks the conda linker); toy runs in 15 s
metadata:
  type: reference
---
2026-10-08, user off-site ("I dont have access to the cluster so you can just test locally"): local builds are fine then.
- Build: ~/SIMPLE_clsrender/build-mac, cmake with ~/.local/bin/devtools-{gfortran,clang,clang++} (micromamba env
  devtools, FFTW there), -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk and SDKROOT set
  to it. The default SDK (MacOSX27.0) fails the C compiler check: libSystem.tbd "unknown architecture" for the conda ld.
  Flags as ~/SIMPLE/build-release (USE_ARCHOPT, USE_OPENMP, -march=armv8.3-a). Full build 3 min 16 s (-j 10, 11 cores).
- Python with numpy+scipy: ~/miniforge3/envs/manifoldem/bin/python3 (system python3 has no numpy).
- Toy: scripts/cls_expansion_proto/fortran_toy/gen_toy.py (fixed seeds = same data as batch2 head002) +
  run_toy_local.sh TOYDIR RUNTAG (scores with score_toy.py and score_frc.py); data in the session scratchpad.
Related: [[cls-expansion-fortran-toy]], [[compile-and-test-means-batch2-fast-gate]].
**Trap (2026-10-13):** the periph8 toy must run with `TOY_APIX=2.135 TOY_MSKDIAM=230` (run_toy_local.sh defaults to 4.27 A /
460). With the defaults the CTF is wrong: class resolution 17 A instead of 4.3 A, which I spent an hour bisecting as a
"master regression" (user: "I think you are chasing ghosts"). Check the import smpd against the reference run first.
