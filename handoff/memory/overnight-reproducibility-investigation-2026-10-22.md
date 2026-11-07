---
name: overnight-reproducibility-investigation-2026-10-22
description: "user 00:20 2026-10-22: 'investigate this overnight' (r1 partition irreproducibility); experiments, job ids, outputs, peer tasks; deliver a morning report"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T04:15:18.409Z
---

The user asked (00:20 2026-10-22): investigate the r1 state-partition irreproducibility overnight. They also said:
- "we shouldn't even be sampling";
- use existing SIMPLE modes, not new strategies;
- the closed-leaning start for the site EM is a concern;
- they approved the third port run on r1.

**Mine** (batch2 `routes/toy55`, outputs `rnasep_klem/repro_r1/`):
- **v2** (`repro_var.sbatch`, `prob_neigh_mode=state` + `prob_athres=0.01`) S/P: 54062656 / 54062660.
- **A, convergence:** `repro_var.sbatch` TAG=v1long, MAXITS=12, greedy; S 54063877, P 54063879.
- **B, how many states:** `repro_k.sbatch` NK=2/3/4/8, seeds = each run's own latent PC1 quantile bins, 6 iterations, greedy.
  - script: 54063881 / 84 / 86 / 88;
  - port1: 54063883 / 85 / 87 / 89.
- **Analysis:** local `scratchpad/viz/cmp2.py` (NMI / ARI / closed-group Jaccard) on the downloaded `*_states.txt`.

**Peers:**
- **Port session:** the third r1 run; merge-up on each r1 run's 8 leaves, then cross-run comparison of the merged partitions (experiment C); greedy projections as the port default.
- **Factorial session:**
  - site-EM pole test (6 starts from latent tails), job 54062884, output `factorial/f692st_r1_poles/groups.txt`;
  - input-pose ranking: done, equidistant from both refined-pose rankings (Spearman 0.92).

**Known before the night:**
- prob_neigh's state labels are greedy already; projection sampling has no effect (v0 vs v1 0.95 from the same start).
- Cross-start: 0.25 / closed-group Jaccard 0.16-0.17.
- `refine=neigh` segfaults in this multi-state focused setup.
- Port pair 0.53 vs script-port 0.25, with no input or field difference found.

**Morning deliverable:** a report doc in `doc/implementation_notes/` plus a chat summary with a recommendation. Related: [[tree-leaves-not-reproducible-across-runs-2026-10-21]].
Progress (00:30): D particle-level width rule = dscore.sbatch (54064637; score_embed per run vs its own 8 stage-B leaf maps -> decisive_merge.py,
report repro_r1/decisive_merge.txt). v2 cross-start 0.24 (no fix); K=2 seeds 0.51 -> stage A 0.56 (ARI 0.67, closed-half Jaccard 0.83) vs
K=8 0.44 -> 0.21; split-then-merge (port, stree15) merges NOTHING in r1 (8 -> 8 every run): the map test cannot reject arbitrary cuts of a
continuum. Pole test: closed found unguided from every low tail. K=2 size coincidence (46129/46425) checked: not eo, row parity or flex A/B.
RESULTS (00:55): K curve (two starts = each run's own latent PC1 bins, greedy): K2 0.51->0.56 (closed J 0.83), K3 0.50->0.51 (0.78),
K4 0.49->0.34, K8 0.44->0.18 (6 it) => r1 supports ~3 reproducible states. Converged K8 (5 it) cross-start 0.24 = different fixed points.
v2 0.24. Particle confusion merge (decisive_merge2.py): paths diverge (script keeps r1_s5 to K=2; ports fold closed into blocks) -> no fix.
Report drafted: doc/implementation_notes/state_tree_reproducibility_overnight_2026_10_22.md (third port run section pending).
LANDSCAPE (01:20, landscape_cmp.py, repro_r1/landscape_cmp.txt): leaf maps matched across runs inside r1's mask vs half-map ceiling.
port1-port2 0.94-0.98 of ceiling (NMI only 0.53!), port3-ports 0.63-0.98, script-ports 0.41-0.89, control r2 0.06-0.31; top-2 map-subspace
cosines >= 0.89 for all four runs, 3rd 0.42-0.96; occupancy shift 6-20 % of spread. => leading 2-D landscape reproducible; finer placement
not. User (01:15): "accept the stochasticity, ensure the landscape is reproducible"; K ladder "not elegant". Overnight tasks sent:
port (landscape tool, final maps, speed), motion (basins on 3 r1 latents), sampling (sampled landscape across 3 runs), factorial (pole set).

**OUTCOME (06:30):** the report `doc/implementation_notes/state_tree_reproducibility_overnight_2026_10_22.md` has a verdict and 7 decisions.
- **r1 node level reproduces.** The integrated port P = 3 runs match within noise at every level, and the delivered closed state is 0.96 / 0.48
  in both. Landscape sampling on the pole axis reproduces, tip included.
- **Merge rule:** see [[pole-merge-line-inside-closed-pair-2026-10-22]]. The recommendation is P = 3 by halves-calibrated pole identity.
- **Depth 1 does NOT reproduce yet.** Children reused the parent's site and their poles degenerate. With their own discovery:
  - the open children have no reproducible site;
  - the closed child finds a reproducible site (Dice 0.90). Its z1 poles match across runs within noise (maps 0.97-0.98) and read
    1.05 / 0.55 vs 0.80 / 0.48, the order/disorder pair;
  - only the argmax seeding over 4 overlapping endpoints is unstable.
- **Next (user's go):** halves-based pole acceptance in site_rank (port est. 4 h), then `site_disc=yes`, then two depth-1 runs.
