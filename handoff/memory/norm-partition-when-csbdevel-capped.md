---
name: norm-partition-when-csbdevel-capped
description: "2026-11-08 USER: \"just use norm\" -> new jobs go to norm (QOS norm, -c <= 36, qsys_partition=norm). History: 2026-10-23 DEFAULT was csbdevel (jobs queue behind our 642 per-user cap); norm also has a 600-core per-user cap reserved for full-data picking/extraction (until ~17:00) -> use -p norm only when told norm is free (<= 36 cpus/job there)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-25T09:52:59.755Z
---

Submit new batch2 work to csbdevel as before, even when it is at our per-user cap: jobs queue behind the 642-core QOS cap.
- norm has its own per-user cap of 600 cores, and the full-data preprocessing (picking, extraction) needs that room. A brief "send new
  jobs to norm" instruction was corrected within the hour.
- Use `-p norm` only when the user or the coordinating session says norm is free. There, keep `-c` <= 36 (36-core nodes) and use
  `qsys_partition=norm` for SIMPLE.
- Never move or edit running or pending jobs or the scripts they use.

- 2026-10-25 go-ahead (coordinator): while csbdevel is capped behind tree workers, small diagnostics (<= 36 cpus: matchers,
  checks, builds) go to norm; tree jobs and their workers stay on csbdevel. Keep ~400 norm cpus free for the user's production
  solve2D (54220038: 48 x 8 workers, ~384 cpus per iteration). Build jobs keep `--constraint=x6342` (norm has x6342 nodes too).

**Why:** 2026-10-23, user via the coordinating session: first "new jobs to norm while csbdevel is capped", then the correction that
norm is capped at 600 and reserved for full-data until ~17:00.
**How to apply:** default csbdevel; switch to norm only on an explicit go-ahead; the 2026-10-25 one covers small diagnostics only.
Related: [[batch2-400-core-limit-and-worker-timeouts]], [[batch2-shared-account-cpu-cap]].
- 2026-10-29 USER: "use norm I have 600 cores there" while csbdevel was capped (406 running / 312 pending). New hetero_tree test jobs
  and their workers went to norm (`-p norm`, qsys_partition=norm, -c <= 36); running/pending csbdevel jobs left untouched.

- 2026-10-29 (user via coordinator: "use norm I have 600 cores there"): csbdevel was at its cap, so NEW jobs went to `-p norm`, including gate trees, with `qsys_partition=norm` in every SIMPLE call so workers land there too; `-c` <= 36. Running/pending csbdevel jobs and their scripts were left alone. sst_gate4 takes QPART=norm for its workers.

- 2026-11-05 USER (IgM thread): "We also have cores on csbdevel, so the investigation needs to be more thorough" -> the earlier
  "IgM stays off csbdevel" rule is lifted for the IgM investigation; audit, replicate and DJ-particle runs use qsys_partition=csbdevel.
  Pass qsys_partition explicitly in every update_project ([[update-project-partition-falls-back-to-bashrc-csbdevel]]).

- 2026-11-06 TRAP: a job's QOS is fixed at submission. `scontrol update Partition=norm` on a pending csbdevel job left QOS=csbdevel, so it stayed blocked on QOSMaxCpuPerUserLimit until `QOS=norm` was also set. Better: submit small IgM diagnostics straight to `-p norm` when csbdevel is at its cap (it was 632/642 running, 37 pending, from other sessions' jobs).
- 2026-11-06 USER (IgM, csbdevel at the 642 cap): "you can shift some jobs elsewhere, norm or unlimited". The `unlimited`
  partition (QOS 432 cpus/user, no time limit, nodes cn084-109/164-167, mostly idle) takes pending SIMPLE workers via
  `scontrol update JobId=<id> Partition=unlimited`; the master does not care where its workers run. New IgM runs may use
  qsys_partition=norm or unlimited.
- 2026-11-08 USER ("umm just use norm lol"), when the coordinator reported 504 cpus on csbdevel near its ~642 cap: NEW jobs go to
  `-p norm`, with QOS norm by submitting there directly, `-c` <= 36 per job, and `qsys_partition=norm` in every SIMPLE call
  so workers land there too. csbdevel stays usable as well; running jobs are left where they are. Relayed to every active
  session.

- 2026-11-08 STANDING (user via the coordinator: "just use norm"): submit NEW jobs straight to `-p norm` (so the QOS is norm),
  `-c` <= 36 per job, `qsys_partition=norm` in every SIMPLE call so tree workers land there too; leave running csbdevel jobs alone.
