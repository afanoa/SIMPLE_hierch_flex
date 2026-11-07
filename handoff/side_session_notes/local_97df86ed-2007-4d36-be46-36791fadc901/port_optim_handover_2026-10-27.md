# hetero_tree / flex_pca native port and performance: hand-over, 2026-10-27

Session: port/perf ("tree-port-optim"), branch `claude/tree-port-optim-de0b3e` at 8c1b66f31, worktree
`/Users/atifao/SIMPLE/.claude/worktrees/tree-port-optim-de0b3e`. Nothing committed, built or run. The integration
worktree `sharp-proskuriakova-b68280` was only read.

Files in this folder:
- `port_keepers_geom_insert_vs_8c1b66f31.diff`: the ported keepers, the same as the branch's working tree.
- `shuffle_optionA_revert20_vs_integration_fa3a5cc7.diff`: option A, one line.
- `shuffle_optionB_sequential_vs_integration_fa3a5cc7.diff`: option B. Both apply to the integration tree's
  `simple_state_mergeup.f90` (SHA256 fa3a5cc7...), not to 8c1b66f31, which has no 50-shuffle rule.

## 1. Inventory of the integration worktree's uncommitted changes

41 tracked files (+6518/-319) and 14 untracked. Every changed line was attributed to the isolated patches in
`scripts/hetero_tree_performance/results/` by exact line matching (script: this session's scratchpad `inv/attrib.py`).
About 3,100 lines match no performance or P-series patch. They are the uncommitted m8c design/experiment layer, row 1.

| # | Change (keys) | Where (lines) | Evidence | Class |
|---|---|---|---|---|
| 1 | m8c design layer: split_ntests Bonferroni re-test (NSHUF_BONF=50), node_motion / strongest_pole / two_poles / motion_gate / motion_inherit / motion_multi / motion_floor, pooled_poles, SIMPLE_SPLIT_CONTRAST / SIMPLE_SPLIT_ROUTE, tree_wiener choice.txt ordering, README | hetero_tree +2466, state_mergeup +383, refine3D commander +189, UI +58, params +18, wiener 3, README 13 | It runs in every m8c build. It causes the 1.81x qtoy slowdown (job 54373629, 20 -> 50 shuffles) | **Not performance work.** Its owner decides whether it is committed (design freeze). Its shuffle cost is item 3 |
| 2 | Projection-geometry reuse, `merge_geom_reuse` | state_mergeup estep, reconstructor tester, params/UI, tree forwarding | Bitwise unit test. 14-group native EM PASS (54372385). Stage timing 54372731: 52.67 / 52.17 / 51.44 s (baseline / off / both), 2 repeats each | **Ported, default off** |
| 3 | Batched multi-target insertion, `merge_insert_reuse`, with its OpenMP and KB-argument compile fixes | reconstructor `insert_planes_targets`, state_mergeup `maps_weighted_shared`, tester, params/UI, tree | Same evidence as row 2. Both on: +2.33% stage speed, no CPU saving | **Ported, default off**, low priority |
| 4 | EM-wrapper audit seam, `merge_insert_audit` / `merge_audit_request` | state_mergeup `em_maps` / `write_em_audit` / `audit_em_wrapper`, refine3D commander +13 | Acceptance tooling for rows 2-3 | **Drop** from production. It stays in the integration tree as test tooling |
| 5 | Fourier-plane storage reuse, `rec_plane_reuse` | image, image_ctf, image tester, matcher_3Drec, UI, tree | OLE 887.0 s vs 697.1 s off (54373389), 27% slower | **Drop** |
| 6 | Packed project subset writer, `node_subset_io` | binoris +104, sp_project(_io), binoris tester, tree +185 | Copy 2.464 -> 0.174 s, but the full chain was 11.726 vs 12.173 s (54357201). Never run natively | **Needs validation.** Port only with the project-lifetime work (section 4b) |
| 7 | Bounded standalone worker pool, `node_worker_pool` / `workers` / `worker_nthr` / `node_worker_mem` | qsys_env / qsys_ctrl / tester, tree +98 | Unit tests only. No native pool started | **Needs validation**, not ported |
| 8 | Independent-state assembly, `state_assembly` / `assembly_mem` / `assembly_state` | volassemble_strategy +477, rec3D_service, flex gateway +6, UI | Modes 0/1 exact. Mode 2 segfaulted (54374185) | **Needs validation**, not ported |
| 8a | `backend = string(queue%get_qsys())` fix and its qsys tester | volassemble_strategy 1 line, qsys tester +26 | Real crash fix for row 8's mode 2. Snapshot 4a29 is unbuilt | **Port together with row 8.** At 8c1b66f31, `get_qsys` has no callers, so this branch has nothing to fix |
| 9 | UI compile/default/placeholder fixes | UI files | Fixes for rows 2-8 | Folded into the ports (both new keys default `no`) |
| 10 | P1 supported extrema, P2 predictive, ancestral init helpers | flex/states, flex/fit new modules, testers, DAG script, highlevel test UI | Research, deprioritised | **Drop** from production paths |
| 11 | P2/P3 callers, full-prior consumer, observation export, population-inference capture, observer v5 | flex_pca basis / embed / fit_driver / fit_types / posterior / probe_fit_update / application / gateway, params +38 | Research | **Drop** from production paths |
| 12 | Stats-writer FLUSH/CLOSE IOSTAT cleanup (writer v7) | `simple_flex_pca_observation_export.f90` (untracked research file only) | Robustness only | **Drop** with row 11. It touches no shared `fclose` |
| 13 | P3 amplitude M-step | probe_fit_estep, 3 lines | | **Parent session owns it.** Not touched |
| 14 | State-map launcher scripts | `scripts/hetero_tree_state_maps.*` (untracked) | 10 offline tests. Native 54358605 / 54361216: payloads identical, 141.6 vs 388.4 s workflow | Tooling. Not ported. Copy into `scripts/` only if wanted |

## 2. Branch state (keepers, uncommitted)

`git diff --stat`: 8 files, +574/-5. The work splits into two reviewable commits:
- **A, `merge_geom_reuse`:** state_mergeup `estep` makes one `new_plane` per particle and band, then `gather` + `put_plane` for
  each live candidate. Also the parameter / parser / validation / UI entries (state_split, state_mergeup, hetero_tree), tree
  forwarding to split and merge-up jobs (subtree jobs inherit the parent command line), and the reconstructor-tester
  `test_projection_geometry_reuse`.
- **B, `merge_insert_reuse`:** `reconstructor%insert_planes_targets` (the final compiled form with both compile fixes), state_mergeup
  `maps_weighted_shared`, and the `em_maps` branch. Also the keys and `test_insertion_geometry_reuse`.
- The ports were adapted to HEAD's `cross_em` / `estep`, which take no design-layer arguments (band fixed at `self%kb`;
  `em_maps` always banded). Arithmetic and order match the single-target path: the same `pf^2 * w` scale, `w * ctfsq`, chunk
  bounds `((c-1)n)/nc+1 .. (cn)/nc`, `compress_exp`, then `sum_reduce` in slot order.
- Two comment lines per commit; no `!<` on the new parameters (the args generator strips them anyway).
- Coordination: the adaptive split/merge session was told by message that these touch `estep`, `cross_em::em_maps` and a new
  `maps` sibling only, not test / decide / eqn_stats / minpop. Its worktree was clean at 8c1b66f31.

## 3. Null shuffles, the 1.81x slowdown: the choice

The cause is design-layer `test()` with `ntests` present, which the tree always passes. All 50 shuffles' maps are built before
the board test. The board threshold uses only the first 20 (max). The other 30 serve only the beyond-EM re-test threshold
`mu + z_crit * sd`, and that runs only when the board test passes (K* >= 2). On qtoy the extra 30 cost 233 s of 519.7 s. Maps,
docking and pair statistics draw no random numbers; only `sshuf` / `rnd_subset` do.

- **Option A: revert to 20.** One line, `NSHUF_BONF = 20`, so the z rule uses 20 shuffles' mean and sd. Its threshold SE is
  0.66 null-sd against 0.41 at 50 (z_crit 3.8). Simulated disagreement with the current 50-shuffle decision, by the
  statistic's distance from the threshold:

  | Distance from threshold (null sd) | 0 | 0.5 | 1.0 | 1.5 | 2.0 |
  |---|---:|---:|---:|---:|---:|
  | z_crit 3.8 | 28% | 18% | 5.8% | 1.2% | 0.1% |
  | z_crit 1.67 | 29% | 7% | 0.3% | 0% | 0% |

  It undoes part of the statistical rule.
- **Option B: sequential, recommended.**
  - All 50 shuffles' subsets are drawn up front in the original order, so the random stream is unchanged; only the 20 the board
    needs get maps.
  - If K* < 2, the other 30 are never built. That part is exact.
  - In the re-test, each merge comparison `best < thr2` builds more shuffles in blocks of
    max(10, ceil(nthr / 2ng)), so every block stays on the per-map path. It stops once
    |best - thr2_n| > 4 * sd_n * spread(n). `spread(n)` is the predicted SD of the 50-shuffle threshold given n:
    sqrt(m/(N n) + z^2 (m-1)(N-2) / (2 (N-1)^2 (n-1))), m = N - n. It matches simulation to within 0-5%.
  - At 50 the band is 0, the same rule as now.
  - Simulated agreement with the 50-shuffle rule:

    | z_crit 3.8 | at the threshold | 1 sd away | 2 sd away | 3 sd away |
    |---|---:|---:|---:|---:|
    | P(decision differs) | 0.24% | 0.01% | 0 | 0 |
    | Mean shuffles | 48 | 38 | 26 | 21 |

    At z_crit 1.67, 2 sd away: mean 20.3 shuffles.
  - Logs gain "from N of 50 shuffles". Board-null statistics are printed over 20 rather than 50. The "shuffle null mean / sd"
    line is dropped.
  - Memory: all subsets are held at once. That is the same as now for ng <= 10 and larger for many groups.
  - It is not literal Besag-Clifford, which counts exceedances: an empirical test at alpha 1/21/m would need about 21m shuffles.
    It is the same idea, stop sampling once the decision is clear, applied to the z rule.
- **Whole-tree effect, either option:** split/merge is 14.5-15.6% of OLE wall and 31.9-38.5% of RNase P, the latter including a
  merge-up already on 20 shuffles. Halving the split step gives at most ~5-15% of whole-tree wall. Not timed.
- **Ownership:** `test()` / `decide()` are the split/merge design area (adaptive split/merge session) and exist only in the
  uncommitted design layer. The patch is drafted, not applied.

## 4. Next targets, scoped

### 4a. Particle reading and preparation (27% of an OLE fit)

- The 27% (59.7 s read + 191.3 s prep of 916.7 s) was measured in a **shared-memory** fit, nparts=1, three passes. Tree fflex
  fits are distributed: fresh workers per part and round, so a RAM cache never outlives a round. Measure the distributed share
  first, on the OLE node fixture with nparts as the tree runs it.
- Prep is likely FFT-bound: three full 720^2 padded FFTs per particle per pass, to use a band of about 22 shells (~1% of the
  plane). That is the real lever.
- Options, from exact to approximate:

  | Option | Size (OLE, 80,599 rows) | Exact? | Saves |
  |---|---|---|---|
  | (i) RAM cache of prepared planes across the passes of one shared-memory fit | 105 GB full, 13 GB per 1/8 partition | Yes, if shift / sigma / half are fixed during the fit (verify) | 2/3 of read+prep, shared-memory runs only |
  | (ii) Band+halo prepared-plane store, written once per node and read by every round's workers | ~2 GB | Yes, values bitwise; storage is a subset of support | Read+prep after the first round |
  | (iii) Compute the band directly, without the full padded FFT | none | No (different arithmetic, like `cache=yes`) | Most of the FFT |

- (ii) needs an audit of every consumer's support first: Cartesian M-step, polar `polar_sample_particle_fused`, PCG, embedding.
  `gen_rec_plane` caps Nyquist only after full-support generation.
- (iii) is like `cache=yes`, which moved 8.5% of OLE assignments while fields moved 0.08%. The OLE partition is that sensitive
  to tiny perturbations, so any non-bitwise preparation change needs an acceptance calibrated to seed-to-seed variability.
  Exact-label gates will fail for reasons unrelated to quality.
- Suggested order: timers -> support audit -> (ii) default off.

### 4b. Project lifetime (site layer)

- `site_layer` defaults to `yes`, which turns `node_proj` off: every node project is root-sized (RNase P: 1,690 files, 64.7 GB).
- Making the site layer row-local means translating at four places in `simple_hetero_tree.f90`:
  - `fs_rank`: `st(self%nrows)` written into the node project -> `st(rows)`. Its site_rank outputs come back local.
  - `fs_climb`: root-sized `start_p*.txt` -> the parent level's local rows.
  - `fs_maps`: global `fs_sets` labels on `stageB/b.simple` -> `st(rows)`.
  - `fs_size`: root-sized rankings.
  - Plus the readers of site results (`site_poles(self%fs_rd(j), self%nrows, ...)`, `fs_output`) map local -> global once.
- Then lift the `node_proj` / `site_layer` guard (`new`, line 171), and port row 6 to remove copy-then-cut.
- **Blocker:** the design layer changes `hetero_tree.f90` by +2466 lines, these routines included. Start only after its commit/drop
  decision, or the work will conflict.

## 5. Validation plan (the user decides when to build and run)

- **Build:** a SLURM job on batch.ncifcrf.gov (never the head node), then `scripts/run_fast_gate.sh build`.
  - `unit_reconstruction` runs both new bitwise tests (C1/C11, banded, 32+3 tail, weights 0 / negative / fractional / >1,
    teams 1 and 4).
  - Plus `unit_ui` and `unit_heterogeneity`.
- **Rows 2-3 on this branch:** the frozen qtoy fixture's recipe passes `split_ntests=132`, which 8c1b66f31 does not know.
  - Use a split fixture without that key (CAHRA 30k C11 stride15, or the qtoy project without it).
  - Run state_split with `merge_geom_reuse` / `merge_insert_reuse` = no / yes / yes+yes, same seed and threads, nthr 1 and 8.
  - Gate: `split_children.txt` identical (`regression_suite.py partition_gate` logic) and the printed EM diagnostics identical.
  - Then one hetero_tree pair (keys off / on) on qtoy, comparing leaf exports with `compare_tree_exports.py`.
  - Timing: palindrome order on one host, as in 54372731.
- **Option B (or A):** apply to a copy of the integration tree, build, run `regression_suite.py --profile exact` on
  `qtoy90_root_split.json`.
  - Gate: labels equal to the 50-shuffle reference in both arms.
  - Also read the native log's "from N of 50 shuffles" and the time against 519.7 / 286.7 s.
  - Add a diagnostic run with split_ntests=1 (z_crit 1.67) and a junk/no-split node, to cover the K* < 2 path.
- **Rows 6-8:** keep their existing unsubmitted harnesses: `native_insert_ablation.sbatch`, the state-assembly protocol JSON, the
  two-part OLE replay for mode 2 with fix 4a29.

## 6. Results, 2026-10-27 afternoon (cluster dir /mnt/beegfs/elmlund/afan/tree_port_optim_20261027)

Builds:
- **int20:** integration tree with the 20-shuffle revert. Job 54383943, fast gate 15/15, simple_exec b9194f57.
- **port:** this branch. Job 54383944, 15/15, bccfc792.
- **intB:** integration tree with option B (sequential null) instead of the revert. Job 54384687, 15/15, 2475bb77.
- **intopt:** int20 plus the reconstruction patch (section 7). Job 54385113.

**Shuffle count, 87 finished tree nodes replayed.** Each node ran its split with two builds back to back in one job, order alternating. The nodes are 25 qtoy plus 32 and 31 OLE base90 nodes (seeds 1 and 2); the 50-shuffle build is the native_callers 385 build.

| | Labels differ from 50 shuffles | Split time vs 50 (OLE, total / per-node median) |
|---|---|---|
| A: 20 shuffles (applied) | 6 nodes: OLE s1 r11121, r11112, r11212; OLE s2 r1122, r1112; qtoy r21112. All 1-5k particles; K* flips both ways | 0.62 / 0.62 |
| B: sequential null | 0 nodes | 0.82 / 0.70 |

- The 50-shuffle build reproduces the original OLE tree labels at every node.
- Exact suite on the qtoy root fixture: PASS for both 20 and 50 shuffles (45 vs 72 s per arm).

**Ported keys (this branch):**
- Labels are identical across off / geom / insert / both at all 56 nodes (8 threads) and on 3 qtoy nodes at 1 thread.
- Time vs off, per-node median at 8 threads: geom 0.90 (OLE) / 0.975 (qtoy); insert 1.01 / 1.12; both 1.06 / 1.15.
- So keep `merge_geom_reuse`; `merge_insert_reuse` is not worth enabling.

**Reconstruction step variants (shared memory, 16 threads):**
- postprocess=no and ml_reg=no save only 0-8% of a step.
- Shared-memory reconstruct3D is not bitwise reproducible run to run: even the `_unfil` halves differ with only postprocess changed.
- The distributed OLE root recB0 spends 162 s in parts and about 400 s in the master. That points at partial-map I/O and 360-box restoration, which cropping removes.

**fflex delivered maps (OLE base90 s1, from file times):** fit 391 s + delivered maps 719 s at the root, 319 + 362 s at r11, 192 + 172 s at r12. The two-group fallback, which reads them, fires at 25 of 34 OLE nodes and 94 of 123 RNase P nodes.

## 7. User direction 2026-10-27: post-flex states at the box the 3D starts at; no redundant reconstructions

Patch `intopt.diff` (in the cluster dir and the scratchpad) against the integration tree:
- New flex_pca key `rec_states_box` (developer; 0 = particle box, default).
  - `deliver_state_maps` reconstructs at that box and registers vol/fsc with the cropped box and smpd.
  - The tree sets it to BOXA (64) for every fflex.
  - The root's first `flex` keeps the native box, because the root focus mask is built from its halves.
- seeds reconstruct at BOXA (they only start stage A); recB at BOXB (88; stage B refs and the children's flex mean).
- postprocess=no on mean, seeds, recB and fs_maps; ml_reg=no plus postprocess=no where only the `_unfil` halves are read (cheap, recB0, refine_rec, pre_rec).
- Deliberately unchanged:
  - recB0 stays native: site discovery uses voxel-count geometry (dilate6 by MDIL voxels), and the children's masks come from it.
  - plain leaf maps (deliverables), mean (the flex working box is unknown to the tree) and stageB (leaf maps).
- Pending: refine3D volrec=no for stA0/stageA, whose maps nobody reads. It skips the per-particle resolution refresh, so it needs a downstream check before use.

## 8. Survey ranking (four read-only surveys, merged)

1. Crop band-limited reconstructions to their consumer's box. This is the largest lever: standalone reconstruction covers about half of tree wall. Done for fflex/seeds/recB in section 7; mean, cheap, recB0 and the site-layer steps need mask resampling first.
2. End phase:
   - submit the merge-up first and run plain/landscape per node as nodes finish (exact);
   - merge-up and split threads should match their 32-CPU allocation (OLE merge-ups ran 8-9 threads for 24-26 min; OLE splits ran 8 threads).
3. Particle prep runs a 720^2 padded FFT for about 22 used shells (exact up to rounding). About 11% of a fit, and it recurs in every reconstruction and refinement.
4. Split / merge-up kernels:
   - longest-first map scheduling (bitwise exact);
   - skip log-only maps under junk_screen=no;
   - full-map reuse in eqn_stats;
   - cross_em total-map hoist;
   - site_em invariant-map and exponential hoists.
5. flex internals: M-step insertion occupancy (about 9 of nthr threads busy), part-reduce streaming, fused embedding Gram, master pre-passes, redundant mean reads.
6. Concurrent site ranking in the default path; shared-memory steps for small nodes; per-step sigma2 reseeding (L).

## 9. Crop step test and map agreement, 2026-10-27 ~12:30

Job 54385342: OLE base90 s1 root, int20 binary, distributed (5 parts), same host for both arms.

| Step | Native box 360 | Cropped | Speed-up |
|---|---:|---:|---:|
| Delivered flex state maps (reconstruct3D m_estimator=flex, 8 states) | 889.7 s | 47.9 s at box 64 | 18.6x |
| recB (2 states), postprocess=no | 255.9 s | 42.8 s at box 88 | 6.0x |

Map agreement (job 54386222). Each cropped map was compared with the native map Fourier-cropped to the same box:
- Box 64 delivered maps: shell correlation to Nyquist (10.7 A) has median 0.999 and minimum 0.975-0.993. Masked CC is 0.990-0.996 for full maps and 0.973-0.977 for `_even_unfil` halves.
- Box 88 recB: shell correlation to 10 A is at least 0.994. Masked CC is 0.998-0.999 for full maps and 0.985 for a half.

Other results:
- Option B exact suite (job 54385252): PASS, 43.6 s per arm.
- The integration tree now carries option B; its state_mergeup sha256 is 34015b0b..., identical to intB.

A/B OLE trees launched 12:08: `ole_s1_base` on intB (2475bb77) is job 54386166, `ole_s1_opt` on intBopt (f48adb40; intB plus intopt.diff) is job 54386167. Both use seed 1, nodepar 8, step jobs pinned to x6342 on csbdevel,norm. Directories: tree_port_optim_20261027/trees/ole_s1_{base,opt}.

## 10. More tree-level patches (2026-10-27 afternoon; scratch git at scratchpad/intdev, files copied to this folder)

Measured from sacct (base90 trees) to choose them:
- Sequential extra-site ranking (`nsite_*`, one job after another until a site splits) costs 33 min on the 692k critical chain and 31 min on OLE's.
- End phase:
  - 692k: plain 4 min, then landscape 3 min, then the merge-up 38 min. The 16-min final_maps gap is already removed in current source.
  - OLE: the factorial site layer runs 1.4 h after all node work. RNase P recipes run site_layer=no.
- OLE merge-ups ran 8-9 threads on 32-CPU jobs. 692k merge-ups ask for 48 threads on a 32-CPU allocation.

Patches (each over the previous; all against the integration tree as of intBopt):
1. `sitebatch.diff` (exact). Once site j gives no split, sites j+1..n are ranked at once. The first of them, in order, that gives at least 2 seed states is chosen, which is the site the sequential chain stops at; otherwise the last site and the two-group fallback.
   - The rankings are independent processes with the same seed.
   - The site layer ranks its own fsite_ dirs (fs_hd), so the extra nsite_ rankings are inert.
   - Restart re-enters from S_SITE for this path (prev_step); batch.txt is written before seed_site.txt.
2. `endphase.diff` (exact). final() submits the merge-up first, then runs plain, landscape and factorial beside it.
   - The merge-up reads only stage-B maps and projects, the leaf table and node masks.
   - final_maps still waits for the plain maps.
   - In subtree jobs, plain and landscape run together.
3. `threads.diff` (numerically close; thread-count chunking).
   - Each step job requests job_cpus_per_task = its nthr.
   - Merge-up nthr is max(nthr, 32, min(48, np/10k)).
   - Split nthr is a thread per 2,500 particles, capped at 32 (OLE 80k: 8 -> 32).

Build intfull = intBopt + 1-3 (job 54387119). Tree `ole_s1_full` launches when it passes. Triple: base (intB), opt (crops only), full (crops + 1-3); seed 1, nodepar 8, x6342.

Not done, with reasons:
- Merge-up concurrent level tests: these change the RNG stream, so decisions can change.
- volrec=no on stA0/stageA: the master reconstruction phase is only 6-22 s of a 200-340 s step.
- Padded-FFT prep: about 4 min per OLE tree, less at box 160, and it is a core kernel.
- Site-layer crops: mask_dice compares masks across sources, and voxel-count geometry would change.

## 11. Bug found by the first tree run (fixed)

ole_nosl_full's root fflex failed: all distributed delivered-map parts stopped with "prg= must come first on the command line".
- Cause: deliver_state_maps set box_crop/smpd_crop on the reconstruct3D cmdline before prg, and part scripts write keys in insertion order. The standalone crop test built its command with prg first, so it passed; the fast gate does not exercise this path.
- Fix: set the crop keys after prg. intopt.diff in this folder is replaced with the fixed version.
- Rebuilds intBopt2 (job 54389317) and intfull2 (54389318). The triple is relaunched as ole_nosl2_{base,opt,full}, site_layer=no, after the first triple (12:28) was cancelled by job ID.

## 12. Timing triple: build hosts (checked 2026-10-27)

All three arms were built by build.sbatch (`--partition=norm --constraint=x6342`, `-DUSE_ARCHOPT=ON`, so `-march=native` targets Ice Lake). None was built on the Skylake login node.

| Arm | Build | Build node | Binary sha (exec = install) |
|---|---|---|---|
| base | intB | fsitgl-hpc119p (cn119) | 2475bb77 |
| opt | intBopt2 | fsitgl-hpc118p (cn118) | d3b69537 |
| full | intfull2 | fsitgl-hpc037p (cn037) | d0bcd766 |

- All three build nodes carry x6342,u467.
- Tree controllers ran on cn119, cn124 and cn132 (all x6342); the workers have SBATCH_CONSTRAINT=x6342.
- So there is no codegen mix between the arms.
- The trees run stage A and stage B with ITA = ITB = 1. The start-point session is moving them to 3+3, so these whole-tree times understate a 3+3 build.

## 13. Tree controller size (recommendation; checked 2026-10-27)

Measured on my 3 OLE controllers (-c 4, 16G), 4 h 21 min in:
- **CPU:** 34-50 CPU-min each, i.e. 3-5 % of the 4-CPU allocation (sstat AveCPU).
- **Memory:** peak 6.0-6.3 GiB.

This agrees with the parent's survey of 35 controllers (8 CPUs: 1-2 % used, MaxRSS 0.6-5.7 GB).

**Recommendation for the production launcher:** 2 CPUs and 16G per tree controller. That frees about 150 CPUs of the 1,242 cap when about 35 trees run.

Not checked: whether any end-phase work runs inside the controller process with OpenMP rather than as a submitted job. If it does, its wall time at 2 CPUs should be watched on the first tree.

## 14. OLE timing triple: results (finished 2026-10-27 22:23)

**Setup:** ole_nosl2_{base,opt,full} on builds intB / intBopt2 / intfull2. Seed 1, nodepar 8, site_layer=no, ITA = ITB = 1, x6342, all launched 13:05:47.

**Whole tree:**

| | base (current) | opt (crops) | full (crops + sitebatch + endphase + threads) |
|---|---|---|---|
| Wall (launcher) | 557 min | 378 min (0.68x) | 474 min (0.85x) |
| Nodes / depth of the longest chain | 31 / 6 | 31 / 5 | 34 / 6 |
| Core-h, all step jobs (sacct, queue-free) | 406 | 277 (0.68x) | 303 (0.75x) |
| Core-h per node | 13.1 | 8.9 (0.68x) | 8.9 (0.68x) |
| Leaves / merged units | 22 / 15 | 26 (1 junk) / 14 | 26 / 14 |

**Why full is slower than opt:** a different tree, not its patches.
- The critical chain is identical through depth 5: 56+38+37+27+87+76 = 321 min (full) vs 57+34+34+36+83+77 = 321 min (opt).
- Full's tree then split once more: node r111111, 94 min.

**End phase** (last node done -> DONE):
- base 25 min: merge-up bound, 21 min.
- opt 47 min: plain 22 min, mostly queue wait, then the 21-min merge-up.
- full 43 min:
  - the merge-up started at once (endphase) and took 9 min (threads), off the path;
  - plain took 28 min (queue-bound);
  - final_maps took 15 min. Its 4 reconstructions (raw, own, tree2, emp2; about 2-2.3 min each) run one after another after plain.

**Per step, (node, step) pairs done in all 3 arms (225 pairs; seconds exit.job -> exit, part-job queue wait included):**

| Step | opt/base | full/base |
|---|---|---|
| fflex | 0.29 | 0.30 |
| recB | 0.39 | 0.39 |
| stA | 0.26 | 0.32 |
| seeds | 0.28 | 0.31 |
| stageA | 0.53 | 0.48 |
| stageB | 0.52 | 0.40 |
| split | 0.96 | 0.77 |
| site | 1.05 | 0.99 |
| flex | 0.96 | 1.02 |
| Total | 0.43 | 0.43 |

- Merge-up: 1275 / 1276 / 550 s.
- "plain" looked 3.5-4.9x slower in opt/full, but it is queue wait only. In node_r222 the master ran 100 s in all 3 arms; the waits were 3+7 s (base, 22:00), 192+244 s (opt, 18:40) and 335+394 s (full, 20:23).

**Queue confound:** median job wait was 120 s (base) vs 21 s (opt) and 14 s (full); p90 about 250-300 s. Step times include part waits, so core-h is the clean compute ratio: 0.68x per node for both cropped arms.

**Absolute walls are high:** nodepar 8 (fast defaults use 16+) and a busy cluster. Use the ratios.

**Leaves do not match across arms:**
- ARI: leaves 0.25-0.27, merged 0.29-0.33.
- Depth-1 (root) ARI is 0.969 base vs opt and base vs full, and 1.000 opt vs full.
- Inside r2, the depth-2 split reproduces in every pair (ARI 0.91-0.97); depth 3 is 0.83 for opt-full.
- Inside r1 (49.8k particles), the split differs in every pair (depth-2 ARI 0.02-0.42):
  - base: site 1 splits 36.9k / 12.9k;
  - opt: site 1 splits 29.9k / 19.9k;
  - full: site 1 gives no split, then batch site 2 splits 30.1k / 19.7k.
- opt and full had identical r1 inputs and still diverged at r1, with no crop difference between them. So the r1 flip is run-to-run variation (non-bitwise steps and a site-order pick at r1, as in the OLE seed study), not a crop effect.
- sitebatch keeps the sequential rule: the first site in order that splits. Its site 1 outcome simply differed.
- Not tested: a same-build repeat (intB twice). It would put a number on the run-to-run floor directly.

**Recommendations from this run (not done):**
- final_maps: run its 4 reconstructions in parallel (about 10 min of end phase).
- plain: submit it per node as nodes finish (the section 8 end-phase item), so it leaves the queue-bound tail.
- Same-build repeat as the leaf-agreement control.
- 3+3 refinement passes (start-point session) will add stage A/B time per node. These 1+1 ratios do not cover that.
