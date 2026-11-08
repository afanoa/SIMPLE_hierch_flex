---
name: sparse-sites-native-port-2026-11-07
description: "flex_pca sparse_sites=yes native port of the sparse-PPCA site pursuit (patch vs e7f7efb80, fast gate 15/15); native is MORE sensitive than the Python prototype (5 % fix 5/5 arms, 5 % null's moving arms found); solvent controls go negative natively, so the per-round check needs the seed-gain floor"
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

**Next:** r131 (both seeds) and 692k closed test on b3 (jobs 55039094, 55039095); then the chain integration.
