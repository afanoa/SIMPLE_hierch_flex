---
name: r131-root-split-flips-with-openblas-variant-2026-11-07
description: TRAP: the r131 gate tree's root split (seed 2, nparts 12) is a near-tie that flips with the OpenBLAS variant linked (libopenblas via -DBLA_VENDOR=OpenBLAS vs libopenblaso by default); every binary reproduces itself byte for byte; build host and code (piece 1) irrelevant
metadata:
  type: project
---

Settled 2026-11-07 on batch2 (math_reworks/gates): mm2 (merge_master2, ee26d9bbf, built on fsitgl-xfer03p without
BLA_VENDOR, links libopenblaso) vs ports/ports3/b0 (ee26d9bbf or ee26d9bbf+ports, x6342 node, -DBLA_VENDOR=OpenBLAS, links
libopenblas): r131 seed-2 partitions agree only .62 after best leaf matching (seed 1: .50), while mm2 == mm2b (rerun) and
ports == ports3 == b0 are byte-identical, and b1 (same as b0 but without BLA_VENDOR) == mm2 byte-identical. The nulls and qtoy2
are identical across all builds. The closed corner is retained in both seeds under either variant.

**Why:** the r131 root split sits on a near-tie; tiny floating-point differences in the BLAS flip it, so gate references are
only comparable between binaries linked against the same OpenBLAS variant.
**How to apply:** compare r131 gate trees only against a reference built with the same BLA_VENDOR choice; record `ldd
simple_exec | grep openblas` in every run manifest; `build_ports3.sbatch` (math_reworks) sets BLA_VENDOR=OpenBLAS, the
merge_master2 recipe does not. Related: [[numerically-close-not-bitwise]], [[hub-simple-infrastructure]].
