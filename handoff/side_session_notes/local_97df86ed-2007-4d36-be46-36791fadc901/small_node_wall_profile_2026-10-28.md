# Where a small node's ~20 min goes: full692 build-7 base tree, 2026-10-28

**Tree:** `/scratch/cluster_scratch/atifao/site_states/full692/base_build_sst7_it3crop_emb_s1`.
- Build 7 = sst7: 3+3 passes, crops, embedding fix.
- The tree started 00:25 and is still running. Nothing in it was touched.

**Sources (read-only):**
- sacct for all 4,878 jobs of the tree (submit, start, end, work directory);
- step markers: `hetero_step.sh` = submit, `.exit.job` = start, `.exit` = end;
- `fsc_*_iterNNN.pdf` mtimes as the pass boundaries;
- the logged "N PART(S) COMPLETED IN x s" lines.

**Scripts:** this session's scratchpad `b7/steps.py` and `b7/nodefacts.sh`.

## Per-step wall (s)

The step job runs from submit to exit, including its own queue wait. Each refinement step runs 2 distributed rounds per pass (prob_tab_neigh, then refine3D); fflex runs 5 rounds (4 flex + 1 reconstruct).

| step (s) | r3222 (6.9k) | r3221 (9.9k) | r3122 (11.3k) | r2321 (12.8k) | r222 (14.7k) | r3121 (21.3k) | r111 (98.8k) | r21 (120.0k) | r (692.4k) |
|---|---|---|---|---|---|---|---|---|---|
| mean | - | - | - | - | - | - | - | - | 104 |
| flex (root) | - | - | - | - | - | - | - | - | 879 |
| cheap (root) | - | - | - | - | - | - | - | - | 144 |
| fflex | 237 | 221 | 243 | 215 | 218 | 263 | 404 | 396 | 839 |
| stA0 | 74 | 85 | 63 | 85 | 86 | 91 | 168 | 177 | 350 |
| recB0 | 46 | 71 | 59 | 45 | 72 | 90 | 134 | 148 | 175 |
| site rank (site) | 26 | 28 | 34 | 34 | 42 | 54 | 173 | 180 | 1220 |
| extra site ranks (nsite_*) | 53 | 119 | 29 | 31 | 34 | 270 | 484 | - | 1144 |
| seeds | 40 | 43 | 44 | 25 | 46 | 53 | 57 | 58 | 136 |
| stageA | 210 | 213 | 198 | 181 | 235 | 240 | 332 | 426 | 815 |
| recB | 41 | 44 | 45 | 45 | 57 | 51 | 61 | 86 | 102 |
| stageB | 228 | 254 | 235 | 255 | 253 | 234 | 542 | 520 | 1350 |
| state_split (incl. null) | 14 | 18 | 19 | 27 | 30 | 40 | 112 | 114 | 303 |
| stageA per pass | 75/66/67 | 77/67/68 | 82/68/48 | 61/53/66 | 92/75/66 | 93/78/66 | 156/74/100 | 159/130/134 | 303/250/253 |
| stageB per pass | 75/77/74 | 84/84/85 | 85/74/74 | 90/73/88 | 75/74/103 | 92/74/66 | 191/172/176 | 194/143/180 | 465/434/442 |
| site ranks run, one after another | 3 | 5 | 2 | 2 | 2 | 7 | 5 | 1 | 2 |
| **sum of step walls** | **969** | **1096** | **969** | **943** | **1073** | **1386** | **2467** | **2105** | **7561** |
| controller hand-offs (step exit -> next submit), total (count) | 198 (10) | 138 (12) | 97 (9) | 109 (9) | 152 (9) | 261 (14) | 301 (12) | 163 (8) | 460 (12) |
| **node span (min)**, first submit -> last exit | **19.4** | **20.6** | **17.8** | **17.5** | **20.4** | **27.4** | **46.1** | **37.8** | **133.7** |

**The same step time, split by cause** (summed over the node's steps, s):

| | r3222 | r3221 | r3122 | r2321 | r222 | r3121 | r111 | r21 | r |
|---|---|---|---|---|---|---|---|---|---|
| step-job queue wait (submit -> start) | 13 | 10 | 6 | 9 | 17 | 13 | 18 | 11 | 11 |
| part-job queue wait | 23 | 27 | 30 | 23 | 28 | 23 | 34 | 26 | 38 |
| part compute (number of rounds) | 151 (22) | 177 (22) | 197 (22) | 200 (22) | 235 (22) | 303 (22) | 795 (22) | 896 (22) | 2290 (29) |
| **master noticing the parts are done** | **487** | **491** | **402** | **360** | **411** | **357** | 371 | 363 | 499 |
| master work between rounds | 117 | 123 | 130 | 130 | 128 | 151 | 264 | 280 | 1457 |
| master or shared-memory work outside rounds (incl. site ranks, split, step start-up) | 178 | 268 | 204 | 221 | 254 | 539 | 985 | 529 | 3266 |

How each column is computed:
- "Master noticing" = the logged phase time minus the round's real span (first part submit -> last part end).
- "Master work between rounds" = the rest of the gap from a round's last part end to the master's next submit.

What the table cannot show:
- The split log has no per-phase timing, so the shuffle null rounds cannot be timed separately. The split total bounds them: 14-40 s on small nodes, 112-303 s on large ones.
- There is no separate "decide" step in this build; K* is decided inside state_split.

## Diagnosis

On small nodes, compute is 151-303 s of an 18-27 min node, and queue wait is only 30-45 s. The dominant cost is the master noticing that a distributed round has finished.

**The evidence:**
- Almost every phase logs "4 PART(S) COMPLETED IN ~30 s", even though its parts ran 4-10 s.
- Across all 857 rounds in this tree, the median gap from a round's end to the master's next action is 25 s.
- The kept beegfs tree (`ole_12707/6_hetero_tree`) logs the same phases at 4.2-10.2 s.

**The cause:**
- `/scratch/cluster_scratch` is NFSv3 with default attribute caching (no actimeo, noac or lookupcache option). The directory attribute cache (30 s minimum) keeps the "file does not exist" answer for the part's sentinel file.
- So qsys_ctrl's 1 s `file_exists` poll (SHORTTIME = 1, the same in `~/SIMPLE_recipe1_src`) sees the sentinel up to 30 s late.
- The tree controller polls `hetero_step.exit` the same way (POLL_SECS = 5). That is the 7-29 s median per step hand-off.

## Levers for small nodes, ranked

### 1. Completion-detection lag on NFS: about 7-10 min per small node (about 40-50 %)

Per node, from 22 distributed rounds, plus hand-offs above about 5 s each:

| Node | Saving |
|---|---|
| r3222 | 10.2 min |
| r3221 | 9.1 min |
| r3122 | 7.2 min |
| r2321 | 6.7 min |
| r222 | 8.3 min |
| r3121 | 8.8 min |

Large nodes save about the same absolute amount (6-8 min), which is a smaller fraction of their time. Options:
- (a) Run tree directories on beegfs. No code change, and the logs above show 4-10 s phases there. Capacity: beegfs is now 62 % used after the deletion, and crops plus map clean-up keep trees small.
- (b) Code: make both polls revalidate the directory before testing the file, for example by opening or listing the directory. NFS close-to-open semantics re-fetch its attributes, which drops the cached "does not exist" answer.
  - It is untested. A 2-minute probe should be run first: a job writes a file on cluster_scratch while another client polls with and without the directory open.
  - Touches qsys_ctrl `update_queue` and the tree controller's step-exit poll.
- (c) Ask the admins for `lookupcache=positive` (or a short `acdirmin`) on /scratch/cluster_scratch.

### 2. Sequential site ranking: 0-5 min per small node, 7.5 min on r111

Small nodes rank up to 7 sites one after another: 26-55 s each plus a hand-off each (r3121: 7 ranks, about 6.5 min). sitebatch.diff in this folder is exact: it ranks the remaining sites at once after the first gives no split. It is not in build 7.

| Node | Saving |
|---|---|
| r3121 | about 5 min |
| r3221 | about 2 min |
| r3222 | about 0.8 min |
| r3122, r2321, r222 | 0 (2 ranks each) |
| r111 | about 7.5 min |
| root | 0 (2 ranks) |

### 3. Refinement passes: 2.5-3 min now per small node for 3+3 -> 2+2, about 1-1.5 min after lever 1

- A small-node pass costs 61-103 s. About 45-55 s of that is the 2 rounds' detection lag, about 6-10 s per round is compute, and about 5-10 s is master work.
- **Risk:** the start-point session's 3+3 evidence comes from a small deep node, RNase P r1211: closed P/R 0.94/0.88 at 2+2 against 0.98/0.98 at 3+3. That is the closed-state hard gate.
- So "fewer passes below N rows" needs a closed-state check before use. After lever 1, the gain halves.

### Not levers on small nodes

- **Shuffle null rounds:** the whole split is 14-40 s, so skipping a second null round saves under 20 s.
- **Queue waits:** step plus part waits are 30-45 s per node at nodepar 4.
- **Step start-up:** about 8-13 s per refinement-type step, about 2 min per node over 11 steps.

## Note on per-level time

- r3222's rows were written at 04:12:52, but the node started at 04:37:33: 25 min waiting for one of the nodepar = 4 slots.
- So a level of n small nodes takes about ceil(n/4) x the node span. Every minute cut from the node span is cut once per wave.

## Lever 1b probe: the fix sees a just-written file at once (2026-10-28, job 54486207)

**Setup:** a 2-node job on norm, with the poller on cn122 and the writer on cn123.
- Directory: /scratch/cluster_scratch/atifao/recipe9/probe.
- The writer wrote a file into each method's own directory 8 times, 8-12 s apart.
- The poller ran `stat` every 0.2 s from before each write, so the "absent" answer was already cached, as in a real master.
- Source: `recipe9/probe/probe.c`.

| method (before each `stat`) | seen | mean lag | max lag |
|---|---|---|---|
| A: nothing (today's code) | 8/8 | 10.55 s | 22.90 s |
| B: open + close the directory | 8/8 | 0.07 s | 0.14 s |
| C: readdir the directory | 8/8 | 0.07 s | 0.14 s |
| D: open + fsync + close the directory (= `simple_sync_dir`) | 8/8 | 0.07 s | 0.13 s |
| E: create + unlink a file in the directory | 8/8 | 0.12 s | 0.18 s |

**The fix used, D:** `call simple_sync_dir(<dir>, status=ios)` before the existing `file_exists` test, in three places:
- `qsys_ctrl%update_queue`, once per poll on `CWD_GLOB`, skipped for qsys_local;
- hetero_tree `step_status`;
- hetero_tree `dir_status`, which also covers `batch_status`.

The diff is `poll_nfs_revalidate.diff` in this folder; recipe9 = recipe7 + sitebatch.diff + poll_nfs_revalidate.diff.

## recipe9 vs recipe7 on qtoy2 (2026-10-28, job 54486894)

**recipe9:** /scratch/cluster_scratch/atifao/recipe9/build_r9.
- Built from recipe7's patched source + sitebatch.diff + poll_nfs_revalidate.diff, on cn128 (x6342), with recipe7's flags.
- Fast gate 15/15.

**Run:** one 2-CPU controller job on norm ran the recipe9 tree, then the recipe7 tree.
- Both used recipe7's prune-test qtoy settings: nodepar 4, maxdepth 1, site_layer=no, oneflex, seed 1, SIMPLE_SPLIT_CONTRAST=maps, RHOMIN 0.03.
- Run directories: recipe9/trees/qtoy_{recipe9,recipe7}_s1.

| | recipe9 | recipe7 (same job, run after) | recipe7, earlier prune test |
|---|---|---|---|
| Wall (tree start -> stop) | **38.7 min** | 56.8 min | 59.9 min |
| Node r (60k rows) span | 22.6 min | 31.3 min | |
| Node r1 (37k) span | 14.3 min | 23.1 min | |
| Node r2 (23k) span | 14.2 min | 22.9 min | |
| Detection lag summed over the node's rounds (r / r1 / r2) | 14 / 14 / 11 s | 394 / 439 / 358 s | |
| Controller hand-off median (r / r1 / r2) | 4 / 3 / 5 s | 22 / 11 / 27 s | |
| stageA passes at r | 45 / 41 / 44 s | 72 / 67 / 67 s | |
| stageB passes at r | 86 / 89 / 83 s | 119 / 116 / 110 s | |
| Part compute at r | 496 s | 503 s | |
| Peak CPUs allocated (+2 controller) | 80 | 80 | |
| Core-h | 15.2 | 17.9 | |
| Leaves | r1_s1 10751, r1_s2 26313, r2_s1 4850, r2_s2 18086 | same | same |

**Timing:**
- Distributed phases complete in 10-44 s on recipe9, against the 30 s / 60 s floors on recipe7.
- Across the tree's 74 rounds, the median gap from round end to the master's next event is 6 s on recipe9 and 25 s on recipe7.
- Both arms have one long r2 hand-off of about 4-7 min. It is not poll-related.

**Leaves are identical:**
- recipe9 vs both recipe7 runs: 0 of 60,000 particles in a different-named leaf, and merged units identical.
- recipe9 and the same-job recipe7 are byte-identical in labels and merged units. The earlier recipe7 run only orders its leaf table differently.

**sitebatch did not trigger:** site 1 gave the split at r, r1 and r2 in every run. Its exactness is therefore not exercised by this test.

**Saving:**
- About 8.7 min per node at 23-60k rows (-28 to -38 %), and wall 0.68x.
- That matches the full692 estimate of 7-10 min per small node.
