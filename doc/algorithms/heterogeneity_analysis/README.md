# Heterogeneity Analysis

Heterogeneity analysis extends the single-map model when one consensus density
cannot explain the particles. SIMPLE provides one continuous model and three
discrete-state workflows, selected by the provenance of the starting model and
particle poses:

- [`flex_pca`](flex_pca.md) estimates continuous variability with a fixed-pose,
  low-rank latent volume model and derives discrete state maps from that latent
  space.
- [`refine3D_states`](refine3d_states.md) refines same-lineage conformational
  states from an existing particle and reference scaffold.
- [`refine3D_states_cavgs`](refine3d_states_cavgs.md) sorts the sub-class
  averages of `cls_expansion` into states, refines them on the class averages
  and hands the labels to particle-level `refine3D_states`.
- [`classify3D_refs`](classify3d_refs.md) classifies particles against a complete
  external reference set, then reconstructs data-derived state maps before
  continuing refinement.
- `hetero_tree` grows a tree of discrete states: every node runs flex_pca, site
  ranking, two refinement stages and a cross-half split test, and its children
  are the groups the test keeps. It has no algorithm chapter yet; its
  distribution options are below.
- `ptcl3D_state_consensus` is a metadata utility in the same category: it
  combines the per-particle state assignments from a file table of projects
  (for example several `refine3D_states` runs) into one consensus `ptcl3D`
  field. It has no algorithm chapter.

## Shared discrete-state model

For `S` discrete conformations, particle `i` carries a state label `s_i` in
addition to its pose, and each state has a map `x_s`. This is the
[refine3D](../refine3d.md) alternation with the reference index extended to
`(state, direction)` pairs: compare the particle with projections of every map,
commit `(s_i, R_i, theta_i, shift_i)`, then [reconstruct](../reconstruction.md)
each map from its members.

Every state is compared with the same whitened Euclidean loss at the same
bandwidth, so the state decision is a likelihood-ratio test:

```text
s_i = argmin_s  min_{R,theta,shift} L_i(s, R, theta, shift).
```

In the probabilistic table the state is chosen by a deterministic argmin over
the heads of the balanced assignment loop
([sampling](../sampling_and_fractional_updates.md)). The projection within
the chosen state is then drawn
stochastically as usual. A shift seed from one state is never reused to rank
another state.

## Shared frequency and coverage schedule

The discrete workflows split the iteration budget into blocks of three and
march the bandwidth linearly in Fourier index:

```text
k_start = max(5, index(lpstart)),      k_stop = min(box/2 - 2, index(lpstop)),
n_blocks = ceil(n_iterations / 3),
k_b = k_start + (b - 1) (k_stop - k_start) / (n_blocks - 2),
lp_b = max(resolution(k_b), lpstop).
```

With defaults `lpstart = 10 A` and `lpstop = 6 A`, the last two blocks run at
`lpstop`. A common band keeps state likelihoods comparable. Without an explicit
`maxits`, the cap is
`clamp(ceil(4 N_active / N_per_iteration), 10, 50)`, approximately four
expected updates per particle.

Fractional updates are balanced over projection-direction bins. A final pass
assigns any active particle with `updatecnt = 0` without changing the last
staged maps. A fresh all-particle reconstruction at native sampling then writes
state maps, half maps, FSCs, and orthogonal reprojections. Single-state
`combine_eo` finalization is not used because merging the half pair would break
the independence required for per-state FSCs.

The common bandwidth and balanced assignment keep the competing maps on equal
statistical footing; the workflow-specific initialization determines whether
that competition begins from a same-lineage scaffold or external references.

## hetero_tree distribution

Every node's steps go through the queue system: `nparts` parts per distributed
step, `nodepar` nodes at a time. With `site_layer=no`, `node_proj` defaults to
yes: each node's projects hold only its own particles, so the controller's step
transitions do not grow with the input (8.3 to 1.7 s per step at a 932k-row
master), and the results stay numerically close to master-sized projects. With
`subtree_max > 0` (opt-in, default 0) a child node of at most that many
particles runs its whole subtree as one job on one node, its steps local. On an
idle queue this costs 2.1-2.5x the worker core-h of the normal tree for no wall
gain (CAHRA 30k: 164 vs 80 core-h, 189 vs 180 min; IgM 100k: 48 vs 19 core-h, 77
vs 76 min), because the subtree job holds its 48 cores through its serial
phases. It pays off only on a contended queue, where one queue wait replaces one
per step; keep it for busy-cluster runs.
