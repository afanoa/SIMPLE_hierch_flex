---
name: confident-cores-replicate-with-the-split-2026-11-09
description: confident closed cores (responsibility >= .9) replicate across seeds exactly where the hard split replicates (r13: cores Jaccard .88, halves .91, share .117 vs .073) and fail where it does not (r131: cores .13, halves .33); selection adds purity, not reproducibility
metadata:
  type: project
---

Three r13-subtree seeds (subtree_repro/r13_s1-3, merge_master2) settled on the test install (prob_beta=1e9 + evidence writer;
node_levers/cores). r13 (same 80,775 rows): settled hard halves 41.8k at closed share .073, Jaccard .908-.912 across seeds
(.83-.84 vs production); cores >= .9 18.8k at .117, recall .25, Jaccard .884-.893, 16,858 shared by all four runs at .1245;
>= .99 5.9k at .160, Jaccard .86. r131 (each seed's own node, root rows): hard halves Jaccard .25-.42, cores .07-.21, shared 158
at >= .9 and 0 at >= .99; one seed's core falls half and half into another seed's two states. r1 (six arm-b settling runs; hits s1-3, s6 with 3 states, misses s4-5): cores >= .9 of the hit seeds 22.6k at .069-.076 (x3.5-3.9),
Jaccard .77-.87 within hits vs .83-.90 for the hard states, .18 hit-miss (a miss seed has no closed-enriched state), 12,308
shared by all six at .0825.

**Why:** the user's rule is "not replicable = redesign"; the coordinator asked whether confident cores reproduce where hard
splits don't. They don't: confidence selection inherits the split's partition.

**How to apply:** deliver "confident cores, with responsibilities" only at levels whose split passes a seed test (r13 yes, r131
no); use the r131 failure as the signal to stop descending, not to select harder. Built on [[responsibility-selection-beats-hard-child-2026-11-09]]
and [[wip-settling-pass-prob-beta-cold]].
