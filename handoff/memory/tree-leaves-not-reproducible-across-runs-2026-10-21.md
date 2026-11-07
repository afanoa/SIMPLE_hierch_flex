---
name: tree-leaves-not-reproducible-across-runs-2026-10-21
description: "TRAP: two runs of the same tree node (f692st r1, 92.5k, script rflex13 vs Fortran hetero_tree) give 8-way partitions that barely agree (NMI 0.24, ARI 0.17); the closed leaf r1_s5 (14.8k, 1.07/0.51) was split over two port leaves (4.9k + 4.9k, each mixed) -> fixed-K leaves of a continuous node are arbitrary cuts"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T02:33:04.135Z
---

2026-10-21: the port's r1 node (`stree_port/tree_f692r1`, hetero_tree_labels.txt) vs the script tree's r1 leaves (`rflex_f692st/final_labels.txt`
1-8), same 92,554 particles and input poses: NMI 0.24, ARI 0.17; no port leaf holds more than 33-57 % of any script leaf. r1_s5's 14,775
particles: 4,915 in port r_s8 (0.90/0.52, plus 3,553 from r1_s6), 4,859 in r_s3 (0.81/0.64, plus 3,763 from the open-ish r1_s3), 2,487 r_s4,
1,618 r_s5. The user found even "slightly less pure" (0.90 vs 1.00) concerning.
LOCALISED (repro_r1.py on existing runs; rnasep_klem/repro_r1/report.txt): whole flex REPRODUCIBLE (script vs port NMI 0.82, closed
group Jaccard 0.77); focused flex not (0.35 same mask; 0.07-0.09 across masks); refinement wanders far from its seeds (fflex -> stage A
NMI 0.12-0.18 in every run; stage A -> B 0.41-0.49) and lands elsewhere each run (final 0.24 r1, 0.12 unit 7; same stage A, other
stage-B mode 0.67). Every run has ONE closed-enriched group (2.3-2.8x, 28-44 % of the node's old closed) but its membership varies
(Jaccard 0.23-0.39; 0.83 when stage A is shared). The split test checks cross-HALF reproducibility of a given partition, never whether
the partition recurs across runs. ROOT CAUSE (repro_e2 + latent_cca.py): refine3D from IDENTICAL fflex seeds+refs reproduces (stage A NMI 0.95, closed group Jaccard
0.98); the flex LATENT is identical across runs with the same mask (canonical corr 0.99-1.00, whole and focused); only flex_pca's K-state
discretisation is random (NMI 0.35-0.61). 8 quantile bins along the shared first latent direction agree at 0.70-0.81, and old-closed
enrichment rises smoothly along it (0.6x -> 1.6-2.0x): a continuum with the closed population at the end. A different mask changes the
latent (canonical corr 0.4-0.7). Deterministic PC1-bin seeds did NOT fix it (seeds NMI 0.44 -> stage A 0.21: refinement amplifies start differences). AUDIT VERDICT
(doc/implementation_notes/state_tree_reproducibility_audit_2026_10_21.md): no bug (identical-start reruns 0.95 / 0.88; no confounding by
direction/defocus/eo/align score, NMI <= 0.02); per-particle evidence is weak: own leaf best 69 %, median margin 1.5-2 nats, only 35 % decided
by >= 3 nats -> the 8-way hard partition is under-determined; width must come from per-particle decisiveness / restart stability.
CONTRAST (factorial session, rnasep_klem/factorial/f692st_r1_sites_reseed): the closed-site 2-class EM ranking in r1 (site mask
r1:states:1f1, r1 stage-B poses) REPRODUCES from 4 starts (r1_s5, port r_s8, port r_s3, seed-free focused flex): Spearman 0.995-0.999,
top-N Jaccard 0.86-0.96 (N 1k-15k), N* 15k, top-2.5k maps 1.16-1.17 / 0.47-0.48; a random start lands on the open pole (sign) -> needs a
closed-leaning seed; the seed-free focused flex is enough. A 1-D site ranking is well posed where the 8-way hard partition is not.
NO CLOSED KNOWLEDGE NEEDED (pole test, factorial/f692st_r1_poles): starts = both 15k tails of focused-latent z1/z2/z3; every LOW tail
-> closed (z1lo, z3lo = the same ranking as the r1_s5-seeded run, Spearman 0.989-0.999, top-2.5k 1.16-1.17 / 0.44-0.48; z2lo a closed
variant), every HIGH tail -> open-type poles (0.39-0.49 / 0.98-1.16). Two-tailed starts on one reproducible direction give both poles;
which pole is the minority/closed one is read from the maps. Input (shared parent) poses: Spearman 0.92 to both refined-pose rankings,
maps ~0.1 closed-arm lower before polishing; port-pose vs script-pose ranking 0.936 (boundary membership follows poses).
Search variants (existing modes) do NOT fix the 8-way cut: greedy projections (prob_athres 0.01) = sampling from the same start (0.95);
prob_neigh_mode=state changes the endpoint (0.50) but cross-start stays 0.24-0.25 for v0/v1/v2; refine=neigh segfaults here.
MASK IS THE LEVER (port, stree_fmask.54067114): port fflex with the SCRIPT's rsmask = script fflex at NMI 0.95 (two seeds 0.93,
closed group J 0.98); with each run's own falloff mask (Dice 0.96-0.98 to the script's) only 0.31-0.36 -> focused flex is ~deterministic
given the mask but hypersensitive to it; run-to-run spread enters via the cheap-children mask. Masks must come from reproducible
quantities (e.g. pole difference maps), judged by downstream sensitivity, never by Dice.
THE FIX THAT WORKS (2026-10-22 ~02:05): K = number of reproducible extrema, seeds = each run's own poles. Pole SET at a shared mask +
shared input poses reproduces across script/port1/port2: closed (Spearman 0.96-1.00, 1.06-1.08/0.47-0.53), open A (1.000), open B
(>=0.989); z2 run-specific. Pole-seeded K=2 stage A (greedy, free assignments) across runs: NMI 0.93-0.95, ARI 0.96-0.98, closed state
J 0.98-0.99 (~36k each) vs 8-way trees 0.18-0.25; maps within noise (1.01 of ceiling), occupancy 0.39/0.61 in every run. K=3:
closed J 0.90-0.96, closed map within noise, open A/B 0.91-1.01, occupancy 0.29/0.34/0.37 in every run. Open: per-run data-derived site masks over-count (broad masks; v2 running);
masks must be reproducible. SETTLED (factorial P1): the site EM is robust to the mask (Dice >= 0.89 -> Spearman >= 0.997); the NODE
DISCOVERY mask (discover_sites2.py, factor split of the refined-state variance, core 0.2) re-derives per run at Dice 0.92-0.94 -> the
same pole. Do NOT use the broad node falloff mask as the site mask (a different, intermediate pole). Axes: anchor them to the poles,
not to eigenvalue order (the top-2 score axes are nearly degenerate and rotate run to run).
**Why:** a state that only appears in some runs is not reproducible; the user's standard is reproducibility.
**How to apply:** never judge a tree change from ONE run against the script tree; compare same-pipeline reruns first (port vs port).
Treat fixed-K leaves of a continuous node as cuts: merge what does not reappear across restarts/halves; the basin (density-peak,
disconnectivity-graph) tree was handed to the motion session for this reason. Related: [[f692st-final-deliverables-2026-10-21]].
