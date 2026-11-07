---
name: gapw-site-acceptance-test-verdict-2026-11-07
description: "VERDICT 2026-11-07 (igm-refs Part 4): GAPw is the per-site discreteness test for sparse-PPCA blocks (9/10 power, 0/10 false accepts on IgM + qtoyflex); per-particle mixture tests only detect non-Gaussianity; fission map tests accept carved continua; GAPw fails on qtoyflex's curved arm continuum (median split DISCRETE .13) and its closed state next to the arm (no call)"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T22:23:53.624Z
---

Oracle session (local_adbcedd9), `hetero_tree_20261028/handoff/igm-refs/VERDICT.md` "Part 4", per-block files `results/part4/`;
GAPw = `fixture_igm/discrete/scripts/dc_gates.py direct` on lever 3's round-1 K = 2 groups of each sparse PPCA block.
- GAPw: IgM arms 9/9 DISCRETE (10 % .03-.24, 5 % .05-.27, 3-state arms split in two included); IgM null blocks .48-1.56 or no
  call; the 5 % non-arm block no call; the null's TRUE angle-median split .47 (not discrete), so the separation is the test's.
- The only calibrated constant is DISCRETE at GAPw < .28 (5 % IgM null median - 3 robust sd), so IgM is partly in-sample; the rest
  are fixed physical / dimensionless constants (5.11 A grid, 15 A lp, 10 A dilation, G2 lobe half-separation >= 10 A, lobe level .5,
  G3 30 deg, G1 z1 >= 3 and r1 >= .5).
- qtoyflex (the RNase P-like control) fails both ways, even with truth groups: closed vs rest = no call (lobes 9 A, G2 fails; the arm
  continuum fills the corridor); the arm continuum split at its true median = DISCRETE .13 (lobes 30 A; a large curved swing leaves
  the chord, GAPw's corridor, empty).
- (a1) held-out 2- vs 1-Gaussian gain: its Gaussian null sits at -.0025 +- .0002 nats/particle, so any non-Gaussian coordinate passes
  (4/10 false accepts). (a2) spikes vs continuum family: misspecified both ways (8/10 false accepts). (c) fission map test: carved
  continuum segments are really different maps (hier-ppca: adjacent qtoyflex arm thirds z 20.5 / 5.2), so it accepts them.

**Why:** the chain sparse PPCA -> one masked round (lever 3) needs a test that keeps continua as one leaf; the seeds cut the null
twin's continuous arms (product .180).
**How to apply:** use GAPw as the site layer's split acceptance on IgM-like short swings, no call = one leaf. Do not trust it on
large curved continua or on a minority state next to a continuum (RNase P's closed corner looks like that): the tree's own splits
must keep carrying the closed state. Related: [[lever3-latent-refs-seeding-verdict-2026-11-07]],
[[sparse-ppca-pursuit-verdict-2026-11-07]], [[discrete-vs-continuous-site-gap-2026-11-05]],
[[cross-half-test-needs-band-above-assignment-2026-11-07]], [[rnasep-closed-state-must-be-retained]].
