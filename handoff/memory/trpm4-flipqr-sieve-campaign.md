---
name: trpm4-flipqr-sieve-campaign
description: 2026-09-18 sieve+abinitio2D runs for trpm4 and flipqr, their working dirs and mask diameters
metadata:
  type: project
---

Queued 2026-09-18 on batch2, partition csbdevel, 128 GB / 10 parts / 8 threads.

Working dirs (both hold a **copy** of the project; originals untouched):
- `/mnt/beegfs/elmlund/afan/trpm4_sieve/trpm4.simple` — copied from
  `/mnt/beegfs/elmlund/processed_test_datasets/trpm4/8_extract/trpm4.simple`
  (read-only, owned by reboulcf). 10488 mics/stks, 1,073,177 ptcls,
  box 180 @ smpd 1.3, cs 2.7, kv 300. Sieve split it into **215 chunks**.
- `/mnt/beegfs/elmlund/afan/flipqr_sieve/flipqr.simple` — copied from
  `/mnt/beegfs/elmlund/afan/flipqr/flipqr.simple`. 2514 mics/stks, 689,169 ptcls.
  Sieve split it into **134 chunks**.

Mask diameters for abinitio2D: **trpm4 = 180 Å**, flipqr = 150 Å. The trpm4 value
had no prior record — measured the automask3D from `trpm4_flex_pcg/automask3D.mrc`
(box 300 @ 1.3), density extent ~191 Å across a 234 Å box, and the user chose 180.
`ncls=300 cache=yes` for both (the user cut 500 -> 300).

Sieve results (2026-09-18): trpm4 **575,795 accepted / 1,073,177** (53.7%, a heavy cut
worth eyeballing in the chunk rejection JPEGs); flipqr **503,881 / 689,169** (73.1%).
Both masters died at the final merge and were recovered — see [[sieve-sigma2-merge-bugs]].

`afterok` is NOT sufficient to chain abinitio2D: when the sieve master fails, the
dependency clears and abinitio2D starts on UNSIEVED data. Add an explicit guard in the
script body that fails closed.

SUPERSEDED 2026-10-02: no SLURM mail at all, see [[no-slurm-email]].

My beegfs area is `/mnt/beegfs/elmlund/afan/` (NOT `/mnt/beegfs/atifao`, which does
not exist; `/mnt/beegfs/elmlund/afan/atifao/` is an old 2025 subdir).

See [[particle-sieving-job-model]], [[flipqr-dataset-setup]],
[[batch2-csbdevel-node-classes]].
