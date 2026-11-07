---
name: rnasep-denoised-stacks-are-phase-flipped
description: "the RNase P denoised stacks (22_denoise_project parts) are ctf=flip in every project that reads them; re-imports must use ctf=flip (klem's --ctf auto already does |CTF|). The broken re-imported maps (contact 0.30) were caused by import_particles resetting state=1 on all rows, not by the flag"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-23T18:19:17.803Z
---

2026-10-16 02:20: re-importing poses into a fresh project (oritab in (stkind, indstk) order + stktab of the ten
denoised parts) paired every image with its pose exactly (checked pixel by pixel) yet reconstructed a wrong map:
same molecule, same frame (corr 0.92 smoothed), arm gone (0.42 = family average), contact 0.30. The family projects
carry `ctf=flip` on these stacks; the import used `ctf=yes`.

Correction 02:45: ctf=flip alone did not fix the maps; the import also resets state (see [[import-particles-resets-state-and-eo]]).

**What the parts are (checked 2026-10-23):** `denoise_raw_particles_part01-10` are denoise_project's REGISTERED phase-flipped RAW
particles, not denoised samples. Each image was rotated and shifted onto its 2D class frame, so every pose refined on them lives in
that registered frame (raw micrograph boxes do not fit them). The interpolation also makes the kx+ky pixel split leak:
[[selected-subset-maps-imprint-inflates-extremes-2026-10-22]].
**How to apply:** `import_particles ... ctf=flip fraca=0.1 kv=300 cs=2.7 smpd=1.464` for the denoised RNase P parts;
check `print_project_field oritype=stk` (ctf= key) of the source project before any re-import; the port's klem
E-step must use |CTF| for these stacks (the prototype's `--ctf auto`). Related: [[rnasep-klem-real-data-setup]],
[[rnasep-692k-project-rows-are-a-permutation-of-stack-order]].
