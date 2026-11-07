---
name: efficient-testing-over-wide-queues-2026-10-27
description: "STANDING (user 2026-10-27): use ALL 1242 cores EFFICIENTLY = keep the cluster full of decisive jobs run in parallel; neither low-value replicate floods nor trimming to near-idle"
metadata:
  type: feedback
---

User, 2026-10-27, a few hours after "test harder / use more cpus": "please prioritize efficient testing (please indicate that to
other chats). We have 1242 cores avil, but we are just waiting in queue becuase of not very useful jobs".

CLARIFIED the same afternoon, after I trimmed every session to ~265 CPUs: "I meant we should be using all the resources of the
cluster efficently." Trimming to near-idle was the wrong reading.

**Why:** after the go-wide push, sessions filled both partitions with replicate matrices (start-point ~700 CPUs, split/merge
~330, my 14-tree matrix queued behind them), so decisive tests waited in queue.
**How to apply:**
- Before each submission, name the decision it informs; seed 1 / one fixture with truth first; second seeds, more fixtures
  and real-data (no-truth) confirmations only after a signal.
- Keep the cluster FULL with decisive jobs run in parallel (independent questions at once), not with replicates whose answer is
  already clear, and never idle while decisive work is pending.
- Relay to every spawned session when it changes. Related: [[autonomous-build-and-submit-2026-10-27]],
  [[batch2-400-core-limit-and-worker-timeouts]], [[price-new-stages-lead-with-lean-scope]].

TRAP (same day): scancel of a hetero_tree controller does NOT cancel its in-flight step jobs; their workers keep being
submitted for a minute. Cancel by job ID, then sweep `squeue -o "%i %Z"` for the tree's run dir and cancel those IDs too.

UPDATE 2026-10-27 ~22:30 (user): "you need to cut these jobs we have a limit of 642 and 600 and I doubt you need to use that
many cores even". The demand (~1000 running + ~1150 pending across sessions, trees at nodepar 16) was itself the problem.
Now: only runs that answer an open decision; demand stays well under the caps; new trees nodepar <= 4 unless a deliverable
needs more; spawned sessions launch nothing without asking the parent. Before adding a run, total `squeue -u atifao` CPUs.
