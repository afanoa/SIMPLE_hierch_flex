# Data-derived node resolution limits (session: tree-data-lp, worktree jovial-jennings-821b3d), 2026-10-27

Cluster dir: /mnt/beegfs/elmlund/afan/tree_data_lp_20261027 (install/, runs/, tools/, offline_lp_54444476.json, score_*.json).
Nothing committed. Total 195 core-h.

## Code (uncommitted, simple_hetero_tree.f90 only)
- Build: 8c1b66f31 + tree_data_lp_on_8c1b66f31.patch; login-node build (fsitgl-xfer03p, Skylake), simple_exec 61807da9.
  Fast gate 15/15, unit_heterogeneity 0.
- The patch also carries ITA = ITB = 3 (the start-point session's 3+3); both arms use it.
- `SIMPLE_TREE_LP=fsc` (env, default off; unset = production behaviour unchanged):
  - stage B lp = node mean's FSC=0.5 resolution (mean/fsc_state01.bin);
  - stage A and flex lp = 1.5 x that;
  - box_crop by refine3D's autoscale(lp * LP2SMPDFAC);
  - oneflex children build a node mean first.
  - The split bands are untouched.

## Verdict: do not adopt; keep the fixed 15/15/10 Å limits; no full trees
- Roots (maxdepth=0, seed 1, off vs on):
  - qtoy2: same within noise (stage-B NMI conformer .426/.429).
  - r1211: closed stage-B recall .937 -> .776 (hard gate fails). Cause: the 2.4 % flex-lp change re-seeded the start-sensitive
    node (seeds ARI -0.02), not a better limit.
  - qtoy null: same verdict in both arms. Committed test 53-54 % (fail); agreed setting rho 0.029 -> K* 1 (pass).
  - OLE null off: committed K* 2, 50.7 % (fail); agreed rho 0.013 -> K* 1 (pass). The on arm died twice on beegfs (below).
- Saved nodes (offline, recB half maps, calibrated to the root means within one shell): in m33 (692k RNase P) the picks coarsen
  with node size (5-10k: 13.8 Å; <5k: 14.6 Å, up to 18 Å; A/flex 18-27 Å). 103 of 122 split-tested nodes fall inside or
  coarser than the 14.7-10.2 Å re-test band, which stage B covers in full today.
- r2121212: off makes one state; on makes a split unrelated to m33's (ARI 0.00) that passes even the rho floor (0.23).

## Traps found
- ~/SIMPLE_smr4's binary predates its SIMPLE_SPLIT_RHOMIN source: it fits contrast but applies no rho floor. Use ~/SIMPLE_smr5o.
- beegfs targets 6202/6203 are full (0 GiB, 0 inodes). Files striped onto them fail, and a rewrite keeps the targets.
  Check with `beegfs-ctl --getentryinfo`; move the step dir aside and rerun.

## 2026-10-29: difference-SNR form (replay on saved half maps; no SIMPLE run)
Dir: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/dsnr_band (dsnr_band.py, list.txt, dsnr_54604768.jsonl; job 54604768,
313 s, 8 CPUs). The definitions are map_merge.py's: envelope-matched S_d / N_d per shell, cross-half. They reproduce its band
values to 3 decimals on the m33 closed pair, the qbp closed arm and OLE r1.
- Inputs: 332 nodes (m33 122, build-7 30 (stage B only), recipe16 r111 12 + qtoy2 21, qtoy2_bp 35, OLE ref6 37, IgM m69M 62,
  8 null roots). recipe14 gates_env and trees9 maps are pruned.
- Where the difference lives (stage-B pairs, classification information M_k SNR^2/(1+SNR)):
  - Real splits: I50 26-34 A, I90 17-25 A. The 20-12 A band holds 3-19 % (m33 6 %); below 12 A ~0 %.
  - Nulls: I90 14.6 A, 20-12 band 32 %, below 12 A 6 %. Nulls sit relatively FINER than real splits.
- Rule (SNR_d >= 1 edge, finest pair, consensus-capped): it holds >= 90 % of the pair information in 279/299 nodes and in
  every named real split. It still coarsens with node size (Spearman -0.7 to -1.0): m33 B 11.7 A at >= 50k -> 23.4 A below 5k.
  It is coarser than the 14.7 A re-test edge at 98/114 m33 nodes.
- Split test: the re-test band (14.7-10.2 A) is where nulls are strong and real splits are weak.

## 2026-10-29: native replays of the difference-SNR band (stage B only), recipe17b
Diff: tree_dsnr_band_on_r17b.diff (this folder; md5 74e05f98234b; simple_hetero_tree.f90 only). SIMPLE_TREE_LP=dsnr: at S_RECB, stage B lp =
finest recB-state pair edge with difference SNR >= 1 (map_merge sums, 3-shell window), capped by consensus FSC 0.143; box by
autoscale. Unset = recipe17b. Build dsnr_band/build_dsnr (login node), fast gate 15/15. The tree's bands equal dsnr_band.py's offline
picks (r111 root 14.6, r111 r1 16.7).
11 roots, maxdepth 0, production keys + split env, seed 1, off vs on; seeds md5 identical in every pair. 304.6 core-h.
- Bands picked: 10.6-21.3 A (qtoy2 root 10.6; r111 root, qtoy null 14.6; m33 r1111221, OLE null, IgM r121 15.6; OLE r1 16.4;
  r111 r1 16.7; r1211, qbp closed arm 18.0; m33 r2121212 21.3). All box 64.
- r1211 closed: state .407/.754 -> .421/.743, group F1 .528 -> .547 (kept).
- qbp closed arm: group F1 .940 -> .938 (kept).
- qtoy null: off passes only via the locality guard (z 5.45, rho .044 > floor); on is refused by the test itself (z 0.98).
- OLE null: rho floor in both arms. qtoy2 root, r111 root, IgM, OLE r1: same K* and partition.
- m33 r1111221 / r2121212: K* 1 in both arms (guard or test); neither reproduces m33's split.
- r111 r1 (user: not worth a split): kept in both arms, weaker in on (z 9.4 -> 5.3).
- Verdict: no gate lost, marginal benefit; off by default under the feature freeze.

## 2026-10-29: ancestral shrinkage (shelved) and the hard band at a deep node (r131)
Diffs (cluster): ancestral/ancestral_shrink_on_r17d.diff (simple_matcher_refvol_utils.f90; SIMPLE_TREE_ANCESTRAL=shrink: per-shell half-map Wiener
shrink of every multi-state matching reference toward the population-weighted consensus); ancestral_dsnr/tree_dsnr_band_on_r17d.diff.
Installs ancestral/build_anc and ancestral_dsnr/build_ad (both gate 15/15; dsnr arm's seeds == off arm's, so off is unchanged).
r131 fixture (44,334, depth-3 context, split_ntests=460, maxdepth 1), 82.6 core-h:
- off: root K* 2 (z 19.9, rho .065); children both split (rho .049 / .039); 4 leaves; corner leaf r2_s1 0.97/0.54.
- dsnr: root K* 2 (rho .091); children at 18.0 / 19.5 A fail z (2.30 / 1.32 vs 3.71); 2 leaves; no corner.
- shrink: root refused by the rho floor alone (z 14.5 passes, rho .0298); 1 leaf; no corner.
- OLE r32 seed pair: deterministic per input (seed 1 == seed 2, ARI 1.0), and recipe17d gives no split even off (z 2.70 vs 3.12).
Verdict: neither band form is safe at deep nodes while the held-out re-test reads 14.6-10.2 A.
