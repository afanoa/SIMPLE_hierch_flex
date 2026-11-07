---
name: tomotwin100-project-setup
description: "How to make the Tomotwin-100 100-species dataset runnable for flex_pca, and the two traps in doing so"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 9671e865-e13f-4298-838f-d873182e5df9
  modified: 2026-08-04T14:17:26.892Z
---

Tomotwin-100 (`~/Downloads/Tomotwin-100.zip`, 6.2 GB -> 7.8 GB extracted) is a **simulated
single-particle** dataset, not tomography despite the name. It is the compositional
counterpart to the conformational IgG set ([[igg-test-project-setup]]).

**Contents:** 100,000 particles = 100 distinct PDB complexes x exactly 1000 each (perfectly
balanced), box 128 at 4.5 A/px (downsampled from 384 at 1.5 A/px), SNR 0.01.
`images/snr0.01/NNN_pdbid_particles.mrcs` (100 stacks), `vols/128_org/NNN_pdbid.mrc` (all
100 GT volumes), `pdbs/`, `init_mask/{init,mask,backproj_snr001}.mrc`.
`gt_latents.pkl` is **(100000,) int64 discrete labels 0-99**, constant in blocks of 1000 —
NOT a continuous latent. Ground truth is therefore purely compositional. Labels also
derivable from `stkind - 1`.

**Import — use the STAR file, not the pickles.** `images/snr0.01/snr0.01.star` carries full
RELION poses (AngleRot/Tilt/Psi, OriginX/Y), CTF, optics and `_rlnRandomSubset`, so
SIMPLE's starproject does the Euler conversion and no hand-conversion of
`combined_poses.pkl` / `combined_ctfs.pkl` is needed.

```
simple_exec prg=new_project projname=tomotwin100
# rewrite _rlnImageName to ABSOLUTE paths first (see trap 1)
simple_exec prg=import_particles projfile=tomotwin100.simple starfile=snr0.01_abs.star \
  smpd=4.5 kv=300 cs=2.7 fraca=0.1 ctf=yes
```

**Trap 1 — relative image paths fail.** `import_particles` runs inside its own
`N_import_particles/` exec directory, so relative `_rlnImageName` entries resolve against
that subdirectory and every stack is pruned ("no stack found for ..."). Rewrite the
`@`-prefixed paths to absolute before importing.

**Trap 2 — `init.mrc` is NOT a data reconstruction.** It is the arithmetic mean of the 100
GT volumes (verified ncc = 1.0000) in simulation units, range 0-21.6. A consensus
reconstructed from the particles ranges 0-0.27. Since flex_pca subtracts `vol1` from the
data, pass a data-derived consensus:
`simple_exec prg=reconstruct3D projfile=... pgrp=c1 mskdiam=450 nthr=8 objfun=cc`
(~1 min for 100k; `objfun=cc` is required — the default euclid path demands a
`sigma2_it_0.star` that does not exist here). It correlates 0.979 with the GT mean.

**Pose validation that should be repeated after any import change:** reconstruct species
000 from its own 1000 particles and correlate against `vols/128_org/000_2cg9.mrc`. Correct
import gives **+0.60** against a cross-species null of **0.15 +/- 0.10**, with unmirrored
beating mirrored (0.60 vs 0.36). A cross-species null this high is expected — all 100 are
blobby complexes in the same box sharing low-frequency mass — so judge by the margin, not
the absolute value.

**neigs must be large.** PCA of the 100 GT volumes inside the shipped mask shows a nearly
flat spectrum: top-20 PCs capture only 41.6% of between-species variance, top-50 75%,
99 needed for 100%. All 100 species are geometrically distinct by k=20 (0 collapsed pairs)
but 2 pairs collapse at k=10 and 10 pairs at k=5. Cache size is the binding constraint,
not compute: `8*(nptcls*ncomp + ncomp^2*nptcls)` bytes = 1.9 GB at neigs=50, 7.5 GB at 100.

See [[flex-pca-embedding-cache-format]].
