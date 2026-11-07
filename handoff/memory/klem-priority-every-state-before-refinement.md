---
name: klem-priority-every-state-before-refinement
description: "User direction 2026-10-15 19:40 on the kernel-latent EM: do not pursue per-state refinement/resolution now; first make sure the latent recovers EVERY state (coverage), judged on toy10 by per-conformer cell ownership and purity (state_coverage.py) and on RNase P by the conformer landscape"
metadata:
  type: feedback
---

User (2026-10-15, after seeing the soft bin maps and asking about filtering): "we dont have to pursue the refinement
now, lets just make sure we can get every state." Earlier the same hour: "the latent seems far more meaningful and is
not a blob so we have to figure out how to push this method so we can reproduce the entire conformational landscape
as the paper" and "so far yes I think this is better than the regular flex".

**Why:** the resolution of the bins is mixture blur + consensus poses (bin 1 vs deposited 18 agrees to 10 A at FSC 0.5
although its even/odd FSC says 6 A); refining spawned states fixes that later. The scientific target is the paper's
55-conformer landscape, i.e. state DISCOVERY and coverage.
**How to apply:** rank experiments by state coverage, not by map resolution; the readout is `state_coverage.py`
(toy55): per conformer the best-cell purity (x base rate), cells owned by majority, recall, and the one-cell-per-
conformer Hungarian purity. Baseline, toy10 full-scale global latent (lat/klem_junk): 10 of 10 conformers own at least
one argmax cell, one-to-one purity mean 0.48 (min 0.36; conformer 1, the central state, is the weakest), best-cell
purity 3-10x base. The full-scale RNase P pose-refinement job (53889289) was cancelled to free cores for the
discovery tests (zoom, sub-class E-step). Related: [[klem-kernel-latent-em-2026-10-15]], [[rnasep-klem-real-data-setup]].

2026-10-15 23:45 (user): "keep working on everything autonomously. The goal is to get as much of the conformational
landscape as possible and all state families. I guess we can't expect significant motion in all those with this
truncated data set." => coverage of states and state families is the objective; per-state refinement and the closed
contact's visibility are parked (15 % of the paper's data); work autonomously, report results as they land.

2026-10-16 morning: see NIGHT_REPORT_2026-10-15.md section "MORNING SUMMARY" (toy55) for the overnight state:
closed state confirmed, node policy settled (classify on the region, align on the whole cell reference), truth-free
region boards recover 8/10 conformers and the docked state, Fortran board + E-step + pose step ported and validated.
