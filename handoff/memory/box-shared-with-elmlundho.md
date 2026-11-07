---
name: box-shared-with-elmlundho
description: "fwl-c142339 is shared: if elmlundho (or anyone else) has jobs running, kill your own jobs and do not queue builds/gates/runs; check `ps -eo user --no-headers | sort | uniq -c` and the load before launching anything"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T00:30:00.000Z
---

2026-09-07: the user stopped a running gate/chain the moment elmlundho's abinitio3D
(nparts=10 + prob_tab workers) appeared on the box: "you should kill your jobs if
elmlundho is running."

**Why:** the box is a shared 80-core machine, not a private cluster; other people's
production runs take priority over my regression gates, builds and timing runs, and my
80-thread jobs plus make -j16 halve their throughput.

**How to apply:** before launching anything on the box, run
`ps -eo user --no-headers | sort | uniq -c` and `cat /proc/loadavg`; if another user has
compute jobs, do not launch, and kill mine if they appear mid-run (anchored patterns,
`pgrep -u atifao -f '^...'`, never a bare name). Queue only when the box is idle, and say
so in the report. Chained overnight scripts must check for other users' jobs before each
stage, not just at launch. Related: [[remote-compute-box]], [[never-signal-simple-processes-by-name]].
