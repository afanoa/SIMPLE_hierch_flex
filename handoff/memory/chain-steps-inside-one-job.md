---
name: chain-steps-inside-one-job
description: multi-step cluster pipelines go INSIDE one sbatch script (switch binaries between steps); never "I'll submit step 2 when step 1 finishes" -- the app can be closed and the chain stalls
metadata:
  type: feedback
---

2026-10-05: the not flex2D control finished 22:44 but its solve3D was only submitted at 10:13 next day, because I had
promised to submit step 2 "once flex2D finishes" with no automatic trigger, and the app was quit overnight. User: "not
sure why flex2D --> solve3D started only 10 mins ago". ~11.5 h lost.
**Why:** my background watchers die with the session; the user also forbids sbatch --dependency chains
([[keep-test-loops-minimal]]).
**How to apply:** write the whole chain as ONE job script whose steps run sequentially (e.g. flex2D on the snapshot binary,
then fix_projfile + solve3D on the latest binary, exporting SIMPLE_PATH/PATH per step, failing closed between steps), or
submit step 2 from inside step 1's script as its last line. Never leave a step pending on my own attention.
Related: [[plain-ssh-sbatch-lacks-gcc-runtime]], [[do-the-requested-job-not-diagnostic-sprawl]].
