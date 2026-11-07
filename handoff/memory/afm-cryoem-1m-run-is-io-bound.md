---
name: afm-cryoem-1m-run-is-io-bound
description: "The 1M-particle afm_cryoem flex_pca run on fwl-c142339 is disk-seek-bound, not compute-bound — threading fixes cannot rescue it without staging the stacks off the spindle"
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T14:42:40.627Z
---

Observed live on 2026-08-12 (`/usr/local/data/atifao/reextract/26_flex_pca`, 3.81M-row project,
~915k state>0 particles, nparts=10 nthr=8, box_crop=64 ncols=128 n_probe_iters=5 box_rec=160):

- All 10 SNR-round workers sat at **~19% of one core each** 43+ min into stage 1 — ~80 cores idle.
- `/usr/local/data` is **one 5.5TB spinning disk** (WD6001FZWX via LVM/xfs), **94% full**; particle
  stacks are **17,668 per-micrograph files** in `108_extract/`. The `rotational storage detected`
  notice caps each process at 2 concurrent reads, but 10 workers x 2 = 20 interleaved streams on one
  head = seek thrash. Box RAM 251 GB (154 GB in page cache) cannot hold the ~390 GB stack set.
- The 1.5 GB project file is re-read by every worker in every qsys round (~10 rounds incl. 5 probe
  iters) — ~150 GB of metadata re-reads through the same spindle.

**Why:** compute fixes (see [[flex-pca-thread-over-particles-fix]]) address the worker CPU envelope;
this run starves before compute matters. Master pre-phase (mean scale etc.) was only ~2 min there.

**How to apply:** stage the live-particle stacks (or a merged stack) onto NVMe — `/home` has 1.5 TB
free on NVMe (nvme0n1/nvme1n1 exist; `/usr/local/data` is NOT on them) — then repoint stk paths.
Failing that: fewer parts x more threads (post-fix nparts=2-3 x nthr=26-40 keeps compute while
cutting concurrent readers), try `SIMPLE_IO_NSTREAMS=1`, and remember every probe iteration is one
more full pass over the stacks (n_probe_iters=5 ~ 5 extra passes).

**RESOLVED 2026-08-12:** rerun as `28_flex_pca` on 23_refine3D_auto's DENSE project (dense rows over
`22_denoise_project/denoise_raw_particles_part*.mrcs`, ~95 GB — fits page cache) with the
[[flex-pca-thread-over-particles-fix]] binary: **932,496 particles end-to-end in ~45 min** at
nparts=10 nthr=8 (was projected 7-10 h sparse). Only the SNR round pays cold I/O (~812 s/worker,
mostly the first 95 GB read); every later round runs cache-warm (probe 194 s/iter, embed 107 s,
states 243 s). Launch env for non-interactive ssh: needs SIMPLE_PATH=~/SIMPLE/build AND
LD_LIBRARY_PATH=/mnt/nasapps/development/modules_elmlund/gcc/16.1.0/lib64 (system libs are too old),
plus mkdir=no. Do NOT just swap vol1 — projfile chooses the stacks.
