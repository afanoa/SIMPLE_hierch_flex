---
name: job-monitor-batch2
description: "user 2026-10-20: 'please have a sophisticated monitor for jobs' -> Monitor-tool job_monitor.py (scratchpad/monitor) + read-only monitor_remote.sh on batch2: one notification per new failure line, abnormal job end, stall, long pending, cap hit, milestone; re-arm every 30 min; TRAP: macOS seq prints big job ids as 5.4981e+07 (use seq -f %.0f)"
metadata:
  type: feedback
---

User asked for it after a silent failure (per-half refinements fell back to consensus poses unnoticed for ~10 min).
**How to apply:** arm whenever long runs are in flight: Monitor(command="python3 -u <scratchpad>/monitor/job_monitor.py",
timeout 1800000) and RE-ARM on every expiry; state persists in monitor/state.json (no repeated events), current table in
monitor/status.txt, planned cancels in monitor/cancelled.txt (add job ids before scancel). Remote side:
routes/toy55/monitor_remote.sh (squeue/sacct/scontrol/grep only - login-node safe); it registers the StdOut of every
running job whose name matches NAMES (extend NAMES in job_monitor.py for new job names) in rnasep_klem/monitor/logs.txt.
Alerts: FAILED|ERROR STOP|Traceback|falling back|no refinement|OOM|empty even/odd|TIMEOUT..., sacct FAILED/OOM/
TIMEOUT/NODE_FAIL/unplanned CANCELLED, running log silent > 45 min, pending > 20 min, partition within 8 of its cap.
Related: [[halves-split-by-eo-break-refinement-trap]], [[klem-board-rewrites-shared-project-trap]].

TRAP (2026-11-06): macOS `seq 54980988 54981039` prints 5.4981e+07 for every id; an sacct-based waiter fed that range waited
on nothing for hours (only the explicitly typed ids resolved). Build id ranges with `seq -f %.0f A B` (or take the ids from
the submit output).
