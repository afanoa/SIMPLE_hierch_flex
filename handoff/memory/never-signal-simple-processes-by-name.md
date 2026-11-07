---
name: never-signal-simple-processes-by-name
description: "pkill by program name kills unrelated production runs (and my own shell); always signal a captured process GROUP instead"
metadata:
  type: feedback
---

`pkill -9 -x simple_exec` / `pkill -f <pattern>` is never safe on a shared compute box.
It has misfired repeatedly:

- **Killed a 100k-particle `abinitio3D` mid-run (2026-08-21).** A *smoke test* of the
  watchdog's stall path executed its `pkill -9 -x simple_exec` cleanup, which is global,
  and took down the unrelated production job at iteration 6.
- **Killed my own SSH shell ~4 times** earlier the same session, because the pattern
  passed to `pkill -f` appeared in the invoking command line itself.

**Why:** signalling by *name* or *cmdline* selects every match on the machine, including
jobs the caller knows nothing about and the caller itself. Correctness here cannot be
recovered by a better pattern -- the set being signalled is simply the wrong set.

**How to apply:**
- In supervisors, put the child in its own process group and signal only that group:
  ```bash
  set -m                                     # bg jobs get their own process group
  "$@" > $LOG 2>&1 < /dev/null &
  JPID=$!; PGID=$(ps -o pgid= -p $JPID | tr -d ' ')
  kill -9 -$PGID                             # ONLY this job
  ```
  This is what `~/watchdog2.sh` on fwl-c142339 does; `~/watchdog.sh` is the old broken one.
- Never smoke-test a destructive code path while a production job is running, even when
  the test looks self-contained -- test the *decision*, not the *cleanup*, or run the test
  on an idle machine.
- If a pattern must be used interactively, bracket one character (`simple_exe[c]`) so the
  invoking command line cannot match itself.

Related instrumentation traps in the same watchdog, fixed at the same time: completion was
detected with `grep -qa "_DONE"`, which matched intermediate markers like `STAGE_A_DONE`
and abandoned supervision of a multi-hour job; and death was inferred from
`pgrep simple_exec` being empty, which is also true in the normal gap between two SIMPLE
invocations. Both are now decided by the child's real exit code via `wait`.
