---
name: sparse-sites-native-port-2026-11-07
description: "flex_pca sparse_sites=yes native port (fork b3ced0021, gate 15/15); native MORE sensitive than Python; seed-gain floor needed; RNase P: 692k block 3 alone = seam-free z_MAP on closed rows (.80), span .85, but the closed leaf pair ranks 75/112 (seam-free 46); real-data stop set by a site-independent floor"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-08T01:05:17.551Z
---

Side session sparse-ppca, 2026-11-07. Port of [[sparse-ppca-pursuit-verdict-2026-11-07]] into flex_pca.

**Where:**
- Handoff: `hetero_tree_20261028/handoff/sparse-ppca/`.
  - `PORT.md`: design, validation table, cost, resume.
  - `sparse_sites.patch`: delivered, md5 837e2639.
  - `sparse_sites_aswritten.patch`: the stop rule as accepted, without the floor.
- Builds: `math_reworks/sparse_port/b3`, delivered (fast gate 15/15, job 55038173); `b2`, as written (15/15).
- Native runs: `math_reworks/sparse_ppca/native/<tag>`, via `scripts/sites_native.sbatch` and `native2npz.py`, which feeds the
  Python scorers.

**Design:**
- Two new modules in `src/main/flex/fit`: `simple_flex_sites_fit` (operator + EM/ADMM) and `simple_flex_sites_pursuit`.
- Data: flex's own whitened observation planes, so ctf=flip gives |CTF| for free.
- Working box `B = 2 round(box smpd / lp) + 4`; constants in its voxels (seed radius 5 for IgM and qtoy alike).
- Only two comment lines: the mandatory `!@descr` headers.

**Traps found:**
- `coords.txt` integrates the 12 nuisance latents under their FITTED prior (`prior_precision`: inv(Gamma_n) on the nuisance
  block). It is NOT the pose profile, though my design note claimed so; site-chain caught it. The site-coupling verdict says
  integrating under the fitted prior invents cross-site couplings. The follow-up `sparse_sites_profiled.patch` adds
  `coords_profiled.txt`: the same solve with a zero nuisance prior block (the Schur complement). The fit stays unchanged.
  COMMITTED 2026-11-08 as fork wip `aac135dfe` (its own build + fast gate 15/15, job 55057869); b4 reproduces b3 to 1e-4,
  692k closed AUROC .802 -> .797. Coupling / merge tests read `coords_profiled.txt`.
- Use the node's `n.simple`. flex refuses its own delivered `fflex/f.simple` (it carries state labels).
- Never nest `contains` inside internal procedures (gfortran).
- A marker test on ssh output catches the login banner: wait on an explicit MARK line instead.

**Findings:**
- Native matches Python to <= .01 where Python found sites.
- Native is MORE sensitive (flex's sigma2 whitening + per-particle normalisation):
  - 5 % fix: 5/5 arms (Python 4/5);
  - 5 % null: the 5 continuously moving arms (R^2 .14-.16 vs T1 .12; disc < 0). Python's "0 sites" was a miss.
- Solvent controls' held-out gains are negative natively (-1.6 to -2.9; Python +0.5 to +4.6), so "2 SE + controls" lets near-dead
  blocks through.
- Fix, with no constant: the grown block must keep at least its seed's held-out gain; drop and continue. Real blocks end at
  1.04-16x their seed, weak ones at 0.07-0.67x.
- Cost: IgM 100k 18-40 min, 13.6 GB (Python 46-125 min, 31 GB).

**Real data (b3, |CTF| via flex's planes; PORT.md "First real-data use", 2026-11-07):**
- r131 root (44k, 3,075 closed): 18 blocks, one site. The best block reads closed at .659 (d' .58); the span .711 only matches the
  seam-free .69-.72.
- 692k root (8,720 closed; 2 h 10 min, 15.6 GB, one node): 16 blocks, one site.
  - Block 3 alone (2 coordinates, 67 A off-centre) reads closed at .802 (d' 1.17), equal to the 30-coordinate seam-free z_MAP
    .798; the span gives .849 (d' 1.43).
  - Tail purity stays low: 5.3 % closed (4.2x the base). The tree's best leaf is 9.3 %.
  - Closed leaf pair `r1111221_s1/_s2` on the seam-free test's 112 pairs (`pair_rank.py`): span d' .48 (seam-free .25), but every
    pair separates more, so it ranks 75 (seam-free 46). No block singles it out.
- **Trap:** the real-data stop is set by a site-independent floor. The tight mask leaves one solvent control, and any local mode
  gains about .004-.008 nats per held-out particle, with no radial trend.
- The strict seed-gain floor dropped nothing on either root (1.5-4.2x the seed), so the strict and 2-SE forms agree there.

- **Trap: the BAND, not the seed size, kills IgM sites at lp 15.** The 10 % fix at lp 15 (box 92) finds no site: seed gain 12.8 vs
  8,172 nats at lp 20.
  - lp 20's 46.7 A seeds at lp 15 still fail (gain 5.5 nats).
  - The same lp 15 lattice and seeds with only the band capped at 20.4 A recover all 5 arms, equal to lp 20 within .007
    (`native/bx10fix20`).
  - The extra 15-20 A shells hold most of a seed's cells and no arm signal, and the shared empirical-Bayes ridge then shrinks every
    seed to zero.
  - User (2026-11-08): "really the band should be data derived". Band scan: box 36 has the gain max at the full band (qtoyflex
    readouts best there), but r131's closed state reads best at 21-26 A while its control floor grows 3x with the band.

**Next:** profiled-coordinate reruns on build b4 (jobs 55047793-5), then the chain integration.
