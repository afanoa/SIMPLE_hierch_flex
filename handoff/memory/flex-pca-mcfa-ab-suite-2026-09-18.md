---
name: flex-pca-mcfa-ab-suite-2026-09-18
description: "2026-09-18 MCFA A/B on the 4-dataset suite: tags mix16 vs mix1 (K=1 = exact plain PPCA) on build ~/SIMPLE_unified_p5 (91c076fc1); final jobs 53236059-66 (see body)"
metadata:
  type: project
---

The unified builds bake the MCFA width at 16 (`COV_EM_MIX`, a7c868d8a removed the env). Commit
91c076fc1 on flex-pca-unified re-adds `SIMPLE_COV_EM_MIX` as an A/B lever only (cov_env_int accepts
>0, so K=1 -- the code's own "reduces EXACTLY to plain PPCA" control -- is the off arm, not 0).
suite_arm.sh takes `EM_MIX=<K>`; suite.sh takes `DEP=<jobid>` (afterok on a build job).

Submitted 2026-09-18 10:05 on batch2: build 53235094 (`~/build_p5.sbatch` -> `~/SIMPLE_unified_p5`,
worktree of 91c076fc1); arms mix16 = 53235095-98 (ribo/splice/trpm4/cnga1), mix1 = 53235099-102.
Runs land in <dataset>/flex_suite/{mix16,mix1}/; score with
`python3 ~/flex_verify/tools/suite_score.py unified6 mix16 mix1` (unified6 = the p4 reference on the
same recipe; expect mix16 == unified6 within the 1e-3 band since the binary only adds an env read).
Reason for the test: prior evidence for MCFA was bundled (10076 AMI 0.180->0.262 came with a delivery
change) or single-row (10028 20k d' 1.46->1.97); the spike ablation showed it neutral on basis
discovery. See [[flex-pca-unified-branch-state]], [[flex-pca-supervised-ceiling-instrument]].

**Resubmissions (2026-09-18):** the first two rounds (53235095-102, then 53235623-28) were cancelled and
their dirs parked as `*_killed_HHMM`: part scripts requested `--cpus-per-task=48` (the master's boosted
nthr leaked into qenv%new; fixed by 4c7f4bb7a on unified / a463b5490 on master-flex-cpus) and the
48-CPU masters plus 48-CPU parts saturated the 642-CPU QoS cap. Final round on the rebuilt p5
(4c7f4bb7a, bin 10:40): masters at `-c 24` (suite.sh), parts at 8 (16 splice):
mix16 = 53235682-85, mix1 = 53235686-89 (ribo/splice/trpm4/cnga1), submitted 10:46.

**Fourth (final) submission 10:57:** the 24-CPU round (53235682-89) was memory-bound (2 TB nodes drained
for maintenance; masters 200-220 GB, cnga1 parts 240 GB, splice parts 125 GB took whole 251 GB nodes;
measured peaks: masters 10-34 GB, parts 8-17 GB). Cancelled, parked as *_killed_1056, resubmitted with
120 GB for masters AND parts (two per node): mix16 = 53236059-62, mix1 = 53236063-66 (ribo/splice/trpm4/cnga1).
suite.sh / suite_arm.sh now carry -c 24 and MEM=120000 everywhere (local copies in scripts/flex_pca_verify
match ~/flex_verify on batch2).
