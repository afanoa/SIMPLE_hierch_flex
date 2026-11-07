# pair_maps FSC-match guard at the hetero_tree gate nodes: stock vs guarded (2026-11-06)

Side session started by the tree-igm-fixture work. Python only: no Fortran change, no commits.

Question: `pair_maps` (src_igm `simple_state_mergeup.f90:5208`) matches group 2 onto group 1 per shell with an unbounded
g = cross / fb. Does guarding it (g = 1 where fb <= 3 / sqrt(N_shell / 2)) change the effect floor (P / N, EFF_MIN 2) or the
locality guard (L / L_ref, TOPQ .05, LOC_MIN .5) at the standing gate nodes?

## Answer

No. The guard changes no verdict, floor or guard, at any gate node. I checked 69 nodes:

- the bc arm (production order: guard, then floor);
- the df replays;
- the base arm;
- the guardfix runs.

It also changes no verdict at the 19 IgM fixture nodes.

- **Why nothing changes.** The guard never triggers at a gate node.
  - The tree's stage-B group maps resolve the whole 10-40 A band. Group 2's half FSC at the band's 10 A end is 0.59-0.94 at the
    bc arm's gate nodes, and still 0.46 for the smallest groups anywhere (2.3k particles).
  - At every band shell of every pair, fb is at least 4.97 x the guard's bar (231 pairs), and |g| is at most 1.15.
  - Stock and guarded therefore agree to the third decimal at every pair.
  - The trap needs a band that runs past the group maps' resolution, as in the IgM prototype's own maps (box 200; fb < .02 at
    11-13 A). The tree's stage-B maps here do not.
- **The gate margins are real, not g noise.**
  - The closed corner holds at both seeds.
    - Seed 2: P/N 3.08 (unclipped 2.99).
    - Seed 1: 2.63 (2.52).
  - qtoy2 r3 and r5 stay merged at 0.22.
  - OLE root pairs are 31-50; depth-1 pairs are 2.57-7.87.
  - The nulls' guard ratios are unchanged (0.204 / 0.234 / 0.237, all merged).
- **The floor does not hold the nulls.** On their root pairs the floor alone would keep the split: P/N 6.1 (qtoynull), 9.9
  (olenull), 40.8 (micnull). Clipped equals unclipped there, so these are real, reproducible nuisance differences. The locality
  guard is the only thing that merges the nulls.

## Gate nodes, bc arm (build_bc: guard then floor; the Fortran logged both values at every node)

P/N columns give the clipped value (the tree's) with the unclipped value in brackets. Guarded columns are shown only where they
differ; they never do. The FSC and g column is per node: the smallest group-2 FSC at the 10 A shell, and the largest |g| in the
band, over its guard and floor pairs.

| node | what it is | children | guard L/L_ref, stock = guarded | floor P/N, stock = guarded | fb at 10 A; max abs g in the band | verdict (stock = guarded) |
|---|---|---|---|---|---|---|
| r131 s2 r1 | closed corner, seed 2 | 12,861 / 12,536 | 1.081 | 3.08 (2.99) | 0.71; 1.03 | split kept |
| r131 s1 r2 | closed corner, seed 1 | 11,344 / 9,541 | 1.190 | 2.63 (2.52) | 0.66; 1.01 | split kept |
| r131 s2 r2 | 0.93 closed-arm substate, seed 2 | 9,850 / 9,087 | 1.129 | **2.08 (2.00)** | 0.60; 1.02 | split kept |
| r131 s1 r1 | sibling of the seed-1 corner | 12,953 / 10,496 | 1.210 | 2.68 (2.55) | 0.68; 1.00 | split kept |
| r131 roots | s2 / s1 | 25,397 / 18,937; 23,449 / 20,885 | 0.922; 0.941 | 7.71 (7.43); 7.97 (7.64) | 0.72 / 0.75; 0.99 | kept |
| qtoy2 r3 | within combo 28 | 6,601 / 3,922 | 1.168 | 0.22 (0.21) | 0.70; 1.00 | floor merges |
| qtoy2 r5 | within combo 101 | 3,562 / 3,045 | 1.142 | 0.22 (0.21) | 0.60; 0.99 | floor merges |
| qtoy2 root | 5 groups | 21,937-6,607 | 0.98-1.11 | 41.9-81.3 (39.5-78.1) | 0.89; 1.02 | kept |
| qtoy2 r1 / r2 / r4 | depth 1 | | 1.04 / 1.01 / 1.25 | 5.17 / 4.27 / 4.74 | 0.59-0.83; 1.02 | kept |
| OLE root | 4 groups | 22,700-17,396 | 1.11-1.25 | 30.9-49.8 (30.3-49.0) | 0.94; 1.05 | kept |
| OLE r1 / r2 / r3 / r4 | depth 1 | r3 11,843 / 7,796 | 1.22 / 1.38 / 1.39 / 0.83 | 7.87 / 5.29 / **2.57 (2.51)** / 6.12 | 0.70-0.83; 1.07 | kept |
| qtoynull root | null | 31,139 / 28,861 | **0.204** | 6.09 (6.08), floor not reached | 0.66; 1.04 | guard merges |
| olenull root | null | 40,610 / 39,989 | **0.234** | 9.90 (9.90), floor not reached | 0.86; 1.01 | guard merges |
| micnull root | null | 32,485 / 24,323 | **0.237** | 40.83 (40.83), floor not reached | 0.83; 1.04 | guard merges |

The replicates give the same picture: every node has stock = guarded and the same verdict. Full tables are in the appendix file.

- **df replays** (lever/splits, build_floor):
  - corners: 3.18 (3.07) at s2, 2.67 (2.55) at s1;
  - qtoy2: r3 0.19, root 41.7-81.6;
  - OLE: depth 1 5.82-6.86, root 29.6-59.0;
  - nulls as in the bc arm.
- **base arm and guardfix runs** (r131_md1_s1 is the same tree as base_r131_s1):
  - r131 s1 r1 / r2: 2.73 / 2.85;
  - qtoy2 depth-2 nodes (guardfix, groups of 2.1k-7.2k particles): P/N 0.26-3.68 and fb at 10 A of 0.46-0.69;
  - olenull root: 0.225 / 12.46.

## Mirror check (stock arm against the tree's own log lines)

- Guard L / L_ref: 111 pairs, median relative difference 0.03 %, max 0.24 %.
- Floor P / N: 55 pairs, median 0.00 %, max 0.21 %. qtoy2 r3 / r5 read 0.218 / 0.220, the same as the tree.
- IgM fixture: guard max 0.54 %, floor max 0.45 %.
- Stage-B labels (ptcl3D states of `split/s.simple`) match the split log's per-state A / B counts at every node.

## Stage-B boxes (read from every header)

- All gate nodes are box 88:
  - RNase P / qtoy2: 2.662 A/px, band = shells 6-23. mskdiam 225 A covers nearly the whole 234 A box: radius 42.4 px + 6 px edge
    against a half box of 44.
  - OLE: 3.903 A/px, band = shells 9-34.
  - micnull: 3.546 A/px.
- The IgM fixture trees are box 164 at 3.988 A/px, band = shells 17-65.

## Sensitivity arms

- **guardf** applies the same 3 sd bar with the soft mask's effective sample fraction, f_eff = mean(m^2)^2 / mean(m^4) = .19-.57
  here.
  - At the gate nodes it is inert as well: the smallest fb / bar is 2.79.
- **g = 1** (no matching at all) changes no gate verdict.
  - It raises P/N by 0.4-25 %, for example micnull 40.8 -> 50.9 and the seed-2 corner 3.08 -> 3.25.
  - On the IgM fixture it flips 3 floor verdicts. fix_root_a groups 1|2 go 2.20 -> 2.63, and the null roots' floor-alone value
    goes 1.94 -> 2.23.
  - So the bounded part of the matching matters near the floor; the unbounded tail does not, at these nodes.

## IgM fixture trees: where the guard touches anything

These nodes are not gates; the check is about mechanism.

- The guard fires at one pair: the real tb100k root's stage-B state pair.
  - At 11.3-12.1 A, fb is 0.010-0.03 and the stock g reaches 2.26.
  - L / L_ref goes 0.480 -> 0.475; the tree logged 0.479. It is merged either way.
- null_c_a2 has 4-5 shells at 16-18 A under the bar. There g is small (0.12), so nothing is amplified.
- Every other IgM pair stays above the bar. The rest-subtracted pairs come within 1.1-1.3 x of it at 15-18 A, and their g stays
  within 0.12-1.10.
- The 39-74 g blow-ups were on the prototype's box-200 maps, not on the tree's own stage-B maps.

## Recommendation for the tree owners

1. **Port the guard into `pair_maps`.** One line: g = 1 where fb <= 3 / sqrt(N_shell / 2).
   - Prefer the effective-sample form, 3 / sqrt(f_eff N_shell / 2), with f_eff taken from pair_maps' own soft sphere. Under that
     sphere, the plain bar is only 1.3-2.3 sd.
   - It is a no-op at all 69 gate nodes and needs no new gate run.
   - It closes the noise path wherever the 10-40 A band runs past the group maps' resolution: small or deep nodes, large particles,
     poorly resolved stage-B maps. That is the IgM prototype's case.
2. **Keep EFF_MIN 2.** With noise-free g, the closed corner clears it by 1.3-1.6x (clipped) and 1.26-1.54x (unclipped), at both
   seeds and in every arm.
   - The tightest kept split is the 0.93 closed-arm substate (r131 s2 r2): 2.08 clipped, 2.00 unclipped.
   - Any change that lowers P/N by more than 4 % merges it.
   - The corner itself is not that close.
3. **Treat the floor as an effect-size floor, not a null control.** The null roots read P/N 6-41; only the locality guard
   (0.20-0.24 against 0.5) stops them. Any change to the guard has to be re-gated on the nulls.
4. **Optional:** log the unclipped sum next to the clipped P/N in the "effect floor:" line.
   - The gap between them is the negative excess that the clipping drops, so it reads the noise in the sum directly.
   - It is 0-8 % of the value at these nodes. On one IgM prototype a1 pair it was 32.2 against 0.66.

## Files (batch2 / beegfs: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/pn_guard_20261106/)

- `scripts/pn_guard.py`: the `pair_maps` mirror. It follows SIMPLE's half-grid FSC / power and `apply_filter`'s shell-0 and
  past-Nyquist rules. It runs the tree's sequence: guard on stage-B state pairs, rank renumbering, then the floor on
  cross-tab-weighted group maps. It has four g arms.
- `scripts/smoke.py` (synthetic check), `scripts/pn_guard.sbatch`, `scripts/nodes.list`, `scripts/igm.list`.
- `out/run1/` (gate nodes) and `out/igm/` (IgM fixture):
  - `nodes.txt`: per node, the box check, labels and verdicts per arm;
  - `pairs.tsv`: per pair and arm;
  - `shells.tsv`: fb, c14, c23, group 1's own half FSC, stock g and both bars at every band shell of every stock pair;
  - `verdicts.tsv`.
- Jobs 54980523 (69 nodes, 7 min) and 54980583 (19 IgM nodes, 11 min): 4 CPUs each, 1.2 core-h in total, MaxRSS 1.3 GB.
- Appendix with every row: `pair_maps_guard_gates_2026_11_06_tables.md`, next to this note.
