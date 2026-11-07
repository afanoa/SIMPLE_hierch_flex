# hetero_tree: coupling the latent model with refinement (verdicts, 2026-11-06)

Status: completed investigation. Nothing from it changes the delivered tree; the one fix it found landed with the master merge
(`ee26d9bbf`). Detailed tables: [flex -> refinement](hetero_tree_flex_to_refinement_2026_11_06.md) and
[M-estimation and ML formulations](hetero_tree_mest_formulations_2026_11_06.md).

## Question

For the 3D tree paper: can flex_pca (the continuous latent model) and the refinement's posterior machinery be integrated more
tightly, in either direction? Constraint set by the user: flex works only at frozen poses, with no latent recalculation or
Woodbury algebra inside refinement.

## Fixtures and test beds

- Gate fixtures: qtoy2 root (truth), the r131 RNase P corner (closed-state reference rows), qtoynull and olenull; the 692k RNase P
  trees m33, base80 and base92 (`rnasep_klem/stree_port/full692_*`).
- Positive control `qtoyflex` (new): qtoy2 with its six arm stops replaced by a 30-map arm continuum (linear blends of the bracketing
  deposited conformers), closed (6 %) and holo discrete. Generator `scripts/cls_expansion_proto/toy55/qtoyflex_maps.py` and
  `run_qtoyflex.sbatch`; data `routes/toy55/simqtoyflex`. The production tree (gate keys, merged build) cuts the arm in two along its
  coordinate and isolates closed poorly: best-leaf F1 0.39 at purity 0.27, holo 0.62 at 0.47, NMI(leaves, kind) 0.255.
- `refine=state_eval` (branch `claude/state-eval`, not on the fork): stage B's own call at the stored poses, each particle's state
  chosen among the states' references at its pose; only state, mi_state and corr are written. It reproduces stage B (qtoy2 99 % of
  labels, NMI 0.656 -> 0.653, no drift over 4 iterations; about a minute per run). An optional per-particle log-prior over the states
  (`SIMPLE_STATE_PRIOR`, `SIMPLE_STATE_PRIOR_ALPHA`) adds a coupling under test. The Python frozen-pose stand-ins drifted below the
  production labels (qtoy2 0.658 -> 0.49-0.52): that decay was their own, not the refinement's.

## Verdicts

| Coupling | Result | Verdict |
| --- | --- | --- |
| Flex prior on the state choice, inside the real refinement (alpha 1, 2) | qtoy2 NMI -0.001 to -0.009; r131 closed enrichment -0.01 to -0.02; toy holo F1 +0.02-0.03, closed unchanged; nulls lose 3-12 % of labels | no |
| Same, calibrated by held-out likelihood (Python) | gain only with site factors; halves on the full node; switches itself off on r131 | no |
| Latent seeding of stage A | seeds start better (qtoy2 0.571 vs 0.309) but stage B ends the same (0.650 vs 0.656) or worse (r131 1.02 vs 1.25) | no |
| State label by evidence summed over poses | r131 closed enrichment 1.450 -> 1.466 (noise), nulls unchanged | no |
| Refit flex's class model from the refinement's posterior | olenull's flex weight goes 0 -> 1 (self-reinforcing) | no |
| Objective reliability weights in flex's fit | about 1 everywhere, anti-correlated with particle SNR | no |
| Flex refit at the refined poses (real flex_pca) | own-latent ceiling 0.477 -> 0.504 on qtoy2, nothing downstream | no |
| M-estimation, soft weights, tree priors, Student-t noise | nothing beats hard labels plus tree_wiener (see the M-estimation record) | no |
| Composed global chart coordinate | seams between children's charts; no minority isolation; merges the 692k closed pair | no (seam-free version under test) |
| Posterior temperature | prob_beta tempers poses only; held-out pose optimum 0.10-0.20, but 0.15 kept the r131 closed corner on 1 of 2 seeds (0.3: 2 of 2, 1.0: 0 of 1) and merged qtoy2's conformer 1/2 pair | keep 0.3 |

## Findings worth keeping

- **Stage B is nearly seed-invariant.** From seeds that agree at NMI 0.03-0.38 it converges to the same partition, so an
  initialisation carries nothing through stages A and B.
- **Per-particle state evidence on real data sits at the nulls' level** (r131 gap 1.8 nats vs 1.4-1.6 on the nulls), and on the toys
  the stage-B mistakes are confident (qtoy2: 32 % wrong, median 0.000 posterior on the true class). The information that separates
  states is in map-level statistics, which the cross-half split tests already use.
- **prob_beta tempers the pose draw only:** `assign_rows` takes the state of the closest reference over all states and draws the
  pose within that state's block. At 1.0 the merged tree lost the r131 closed corner.
- **The split step is not reproducible on a continuum:** identical stage-B input gave different children on qtoyflex (holo F1 0.62
  vs 0.49).
- **The r131 closed corner sits at the readout threshold** under every setting tried (closed arm 0.91-1.03).
- **Site factors are the only refinement -> flex signal:** they lift the latent's view of the states (qtoy2 truth ceiling NMI
  0.48 -> 0.62, r131 closed AUROC 0.65 -> 0.69).

## Where things are

- Python prototypes: `scripts/hetero_tree_proto_latent/` (untracked, main checkout), on batch2 under
  `hetero_tree_20261028/math_reworks/latent/` (runs, dumps, caches, `seval_*`, `global/`).
- Gate arms: `math_reworks/gates/{mm2,mb1,pb15}_*`, `mm2_r131_s1`, `pb15_r131_s1`, `mm2_qtoyflex`; readouts `gates/<arm>_readout.txt`.
- Traps met: flex_pca refuses a project carrying several state labels (make one population with
  `simple_exec prg=selection oritype=ptcl3D infile=<0/1 per row> mkdir=no`); `selection` refuses node-local projects whose stack
  segment is the parent's; `simulate_particles` refuses a class with an odd particle count.
