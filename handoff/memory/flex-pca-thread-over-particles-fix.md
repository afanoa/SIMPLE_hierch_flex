---
name: flex-pca-thread-over-particles-fix
description: "SNR/columns/mean-scale stages rethreaded over particles on 2026-08-12; columns is provably bitwise, and SIMPLE_COV_DUMP_ACC is the A/B instrument"
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-12T19:29:40.835Z
---

On 2026-08-12 the three flex_pca stages that serialized per particle were rethreaded in
`src/main/flex/simple_flex_pca_columns.f90` (uncommitted, like the rest of flex):

- `estimate_snr_volume`: parallel-over-particles with per-thread scratch reconstructors and
  per-thread var/dens accumulators, folded in fixed thread order. Reproducible at fixed nthr;
  changing nthr moves sums by ~1e-5 rel (float32 reorder), which feeds only the selection argmax.
- `accumulate_covariance_columns`: two-phase blocks of nthr particles — phase A parallel-over-particles
  builds caches/gathers into per-slot arrays; phase B parallel-over-COLUMNS scatters in particle order.
  **Bitwise identical at any thread count by construction** (each column owns a disjoint slice, same
  per-voxel addition order). `backproject_columns` was deleted; its body lives in the contained
  `scatter_one_particle`.
- `estimate_mean_scale`: per-thread partials, ordered fold. dp drift vanished in the float32 filt cast
  on IgG — whole pipeline stayed byte-identical, but that is dataset luck, not a guarantee.

Validated on IgG 10k (box_crop=64, ncols=24): cols accumulator dumps byte-identical across
base/new and nthr=8/3; column selection stable; SNR drift max_rel 1.3e-5. Laptop nthr=8 timings:
SNR 19.1→11.6 s, columns 23.2→11.6 s.

**How to A/B any future rethreading:** `SIMPLE_COV_DUMP_ACC=1` makes shmem runs write
`flex_pca_snr_dump.bin`, `flex_pca_cols_dump.bin` (raw accumulators, cmp-able as float32 streams
after the 36-byte header) and `flex_pca_columns.bin`. Run with `column_sampling=lowfreq` to test the
columns stage in isolation from SNR rounding. Kill the run once the cols dump appears (it is
tmp+rename atomic). A/B binaries must keep basename `simple_exec` — argv[0] is validated.

Deploy note: the fix exists only on the Mac worktree; fwl-c142339 runs its own tree at
~/SIMPLE and needs the file synced + rebuilt + `cp build/production/* build/bin/` (see
[[flex-pca-ab-testing-traps]] trap 1). Related: [[afm-cryoem-1m-run-is-io-bound]],
[[flex-pca-dtilde-and-speed]].
