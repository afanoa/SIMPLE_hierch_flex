---
name: batch2-400-core-limit-and-worker-timeouts
description: "user 2026-10-21: NO 400 limit - >= 642 csbdevel + 200 norm; 2026-10-27: 1242 = 642 csbdevel + 600 norm shared by all sessions, under-use is the problem; starved klem boards time out waiting for workers and a failed region stalls its family job for 12 h"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T17:57:37.856Z
---

2026-10-19: "we need quicker tests to establish viability. I also need to process other data so we should
avoid using more than 400 cores for now." Six half trees had filled the 642-core QoS cap (288 cores
were six idle 48-core chain masters).

**Why:** the cap is shared with the user's own processing; and a klem board whose SLURM workers sit
pending dies with "timed out waiting for worker N reply" (simple_klem_distr), the region writes no DONE
marker and f692_regions_par3.sh waits 720 min for it, so starving my own jobs stalls the chain.

**How to apply:** count cores before every submission (`squeue -u atifao -h -o %C | awk '{s+=$1}'`),
keep the total under 400, size chain masters with MC=16 (synth_tree_prep_full/split_halves), run
small pairs (NP=5, PARMAX=2) and skip level-2 zooms (NL2=0) when only level-1 leaves are read
(merge_cross_tree.py and consolidate_states.py use level1 only). A failed region is resumable:
move its board dir aside and resubmit f692_region_one3.sbatch with the family's env; the family's
wait loop picks up the DONE marker. See [[batch2-shared-account-cpu-cap]],
[[no-4m-synthetic-tests-too-expensive]], [[klem-cross-tree-merge-design-2026-10-19]].

UPDATE 2026-10-19 ~14:40 (user): "we can use up to 642" - the 400-core limit is lifted; the QOS cap (642) is the limit.
Other sessions' klem_dtree/klem_dtboard jobs can fill it within seconds (416 cores at 14:43), so submit board-heavy
chains with LOCALDEV=1 (family = one 48-core job, no SLURM board workers that time out while pending).
TRAP: never rm a run directory right after resubmitting the job that uses it - the job can start between squeue and rm.

2026-10-19 later: the user's own dataset processing (ribosome_10509 preprocess, 20 workers) takes priority over the
400-core guideline ("the ribosome processing takes priority"); the guideline caps experiments, not real processing.
Still tell the user the account total when processing pushes it past 400.

2026-10-20 (sacctmgr show qos): per-user CPU limits are per QoS and partitions have their own QoS: csbdevel 642
(group partition, 2,928 CPUs), norm 600 (default partition, 4,800 CPUs / 101 nodes, DefMemPerCPU 20 GB -> always pass
--mem), unlimited 432 (no walltime cap), short 170 (30 min); the user association has no extra cap -> csbdevel + norm
can run side by side (~1,240 cores). SIMPLE workers follow SIMPLE_QSYS_PARTITION / qsys_partition, not #SBATCH -p.
When one QoS is full (Reason QOSMaxCpuPerUserLimit) and the other has room, PENDING jobs of mine can be moved: `scontrol update JobId=<id> Partition=csbdevel` (or norm); 8 synthetic klem workers moved norm -> csbdevel 2026-10-20 11:38 and started at once.
The 'unlimited' partition (QoS 'infinite', 432) is NOT usable as extra capacity: 2026-10-20 jobs submitted with -p unlimited (header or command line) came back as Partition=norm, QOS=norm (rerouted, probably because the walltime fits norm) -> the usable per-user budget is csbdevel 642 + norm 600.
TRAP 2026-10-20 15:50: `scontrol update JobId=X Partition=norm` alone keeps the job QOS=csbdevel -> it stays pending (QOSMaxCpuPerUserLimit) under the csbdevel per-user cap; also set `scontrol update JobId=X QOS=norm` (the swing sim started at once).

2026-10-21: "why are so many of my jobs on norm. The queues are longer there generally" (after an earlier "you can use norm if
needed" that I had turned into a default). Submit to csbdevel (and QPART=csbdevel for distributed SIMPLE workers); norm only
when csbdevel is at the cap, and say so when I do.

UPDATE 2026-10-21 14:30 (user): "forget about the 400 limit. You can use at least 642 on csbdevel and 200 on norm if you want" - the 400-core guideline is gone; budget = csbdevel QOS 642 + 200 on norm (submit with -p norm and --mem; a moved pending job also needs QOS=norm). A peer session relaying "under about 400" is superseded by this.

2026-10-21 later: "core budget is 642 csbdevel + 200 on norm" -> the 400-core cap above is superseded; plan jobs within
642 csbdevel + 200 norm across ALL my sessions (count other sessions' jobs before scaling up).

2026-10-23 (user, relayed): while csbdevel is at the per-user cap (642 running, jobs pending on QOSMaxCpuPerUserLimit), send NEW jobs
to norm (`sbatch -p norm`, `qsys_partition=norm` for SIMPLE trees/runs); leave running/pending jobs alone and copy any script a
running job uses rather than editing it. norm nodes have 36 cores (csbdevel 48), so keep per-job cpus <= 36 on norm; both allow 5
days and >= 256 GB per node. Full-data preprocessing keeps priority on norm; report it if my norm jobs pend on QOSMaxCpuPerUserLimit.
CORRECTED the same day: norm has its own 600 per-user cap and full-data picking/extraction needed it (~502 in use, 14:40-17:00),
so new jobs went back to csbdevel, queueing behind the 642 cap. Ask before using norm; the coordinator says when it frees up.
CORRECTION same day: norm ALSO has a per-user cap (600 cores; QoS norm). Full-data preprocessing (~464 cores) needs that room, so new
jobs go to csbdevel and queue behind its 642 cap until norm frees (relayed: until ~17:00). Already-running norm jobs stay. Before any
norm submit, count the account's norm cores against 600, not just the idle nodes.

UPDATE 2026-10-23 (about 23:00, user): "you can use 600 + 642 cores as you want". The full per-user budget, csbdevel 642 + norm 600 =
1,242 cores, is available for the overnight / autonomous work across all sessions. Still count our cores before submitting, and keep
the full-data processing's share (e.g. its 400-core norm sieve) running.

2026-10-24 about 21:30, user: "we have a lot of cpus available. please proceed autonomously". So for directions the user has approved,
don't hold each step for a price approval:
- set -t to about 2x the expected time, and report costs as they land;
- run independent checks in parallel, not one after another;
- still stop by id and re-price anything heading well past its estimate ([[price-new-stages-lead-with-lean-scope]]).
- norm also carries the user's full-data 2D (384 cores) right now, so the default stays csbdevel ([[norm-partition-when-csbdevel-capped]]).

2026-10-27 (user): "we can use 642 cores on csbdevel for development". This is ONE pool shared by every concurrent
session (parent amplitude A/B, adaptive-split-merge, tree-port-optim). Rule sent to both: before each sbatch sum the
account's running+pending csbdevel CPUs (`squeue -u atifao -p csbdevel -h -o %C`, workers included) and submit only if
total + request <= 642; single submissions <= ~320 unless the queue is empty; small diagnostics (<= 36) may use norm.
SUPERSEDED same day (user, 2026-10-27): "we need to use more cpus for testing ... We have acces to 1242 cores yet I have only
alloated 400." Budget = 642 csbdevel + 600 norm, shared by all sessions; under-use is the problem. No per-submission cap; only
the per-partition running+pending check. norm is OK for tests incl. workers, but pin --constraint=x6342 (+ SBATCH_CONSTRAINT)
because builds are -march=native on x6342 and norm also has x6150 nodes.
PACING RULE (2026-10-27 afternoon, sent to all sessions): do NOT count other sessions' PENDING jobs (an unrelated simple_dec
chain had 360 CPUs pending on norm and made every pacer read "over budget" while cores sat idle). Submit while
(account RUNNING CPUs + own pending + request) <= 642 csbdevel / 600 norm; SLURM enforces the account caps and queues the rest.
CAPS ARE QOS PER USER (checked 2026-10-27): csbdevel MaxTRESPerUser cpu=642, norm cpu=600; atifao's association allows only
QOS `normal` (the lab's norm_ncifelmlucryosparc cpu=1500 belongs to the cryoSPARC account). csbdevel HARDWARE had ~1,470 idle
CPUs while 233 of our jobs waited on QOSMaxCpuPerUserLimit -> more capacity needs an admin cap raise / QOS grant.
scontrol top is NOT permitted; scontrol hold/release of own pending jobs works.
TREE CONTROLLERS ARE ~1-2 % BUSY: 35 controllers held 190 CPUs (8 each) at 1-2 % use (41 CPU-min / 5.5 h). Use 2 CPUs + 16G
for hetero_tree controllers (MaxRSS 0.6-5.7 GB); workers keep nthr via step args.
