---
name: cls-expansion-stage1-ready-2026-10-12
description: "cls_expansion stage 1 (conditional averages, no solvent prior) is on origin/master as c8662fb2a (2026-10-13); stage 2 = solvent prior from cls-expansion-solvent"
metadata:
  node_type: memory
  type: project
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-13T02:10:32.062Z
---

Stage 1 of the two-stage cls_expansion push: worktree ~/SIMPLE_stage1, branch cls-expansion-stage1, committed
2026-10-13 as c8662fb2a rebased on origin/master eed4f9609 (11 files, +865/-956); batch2 build + fast gate 15/15
passed; PUSHED to origin/master 2026-10-13 (fast-forward eed4f9609..c8662fb2a) on the user's "push it". Conditional averages P + U z_s, even/odd from independent
half models, FRC filter, residual junk axis kept, rank 16, data-derived band; soft weights/cross-fit/refine_inpl/restore
removed; cls_expansion_frcs.bin registered as frc2D (fixes the stale per-parent frc2D that solve3D_cavgs reads for its
lp ladder - a master bug too). Verified on the Mac only: unit suite, 3-class toy matches master's contract, 1 vs 2 parts
bit-identical. Real-data runs 2026-10-13 in /mnt/beegfs/elmlund/afan/clsexp_stage1/{rnasep,10180,not}
(5 classes each, ncls=100). Open: the stream (simple_stream_stage_solve3D start_clsexp) passes ncls per class
but cls_expansion reads the total - stream-side fix, flagged to the user, not ours.

Stage 2 (2026-10-13): 2.1 circle mask PUSHED as 8e727d8ae (on c6e344c6b); 2.2 solvent prior (Otsu core + spread,
lambda 1 constant, no parameter) reworked onto SIMPLE's pcg_solve / image FFT / otsu / image_bin grow_bins+cos_edge
(cos edge must be CENTRED: grow rdil/2, falloff rdil - else support too wide, head purity 0.91 -> 0.81; PCG rtol/xtol
1e-7) plus one ctf_mean/shell_frc_dp/spd_inverse; rebased as 498ff04e7 (2026-10-14).
**2.2 DROPPED by the user (2026-10-14)** after the kernel-route test (subcls_avg=kernel, master 6d56e6da1): the prior was
neutral on real kernel averages (RNase P / 10180 / NOT montage near-identical). Not pushed; branch cls-expansion-mask
keeps it locally. automask2D as the core lost the beyond-mskdiam arm on the realistic toy (0.95 -> 0.75) - rejected.
Realistic toy needs TOY_APIX=2.135 TOY_MSKDIAM=230.

**Why:** the user wants the push in two stages and to see real data first.
**How to apply:** commit only when asked (no co-author trailer), never push without explicit permission; stage 2 =
the solvent prior from branch cls-expansion-solvent (~/SIMPLE_clsclean). Related: [[cls-expansion-solvent-prior-2026-10-09]].
