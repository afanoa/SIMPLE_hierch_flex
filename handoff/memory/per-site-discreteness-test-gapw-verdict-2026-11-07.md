---
name: per-site-discreteness-test-gapw-verdict-2026-11-07
description: "VERDICT 2026-11-07 (Part 4 of igm-refs): GAPw on round-1 groups is the per-site discreteness test (IgM 9/9 arms, 0/9 null/non-arm) but fails qtoyflex both ways (closed missed with true labels; curved arm continuum split at its true median called DISCRETE .13); per-particle 2-vs-1 / spikes-vs-continuum tests fail (non-Gaussianity detectors at N 1e5); fission map test accepts carved continua"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-08T00:19:54.495Z
---

Handoff `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-refs/VERDICT.md` section Part 4; per-block files results/part4/.
Scripts `disc_a.py` (+ sbatch / list), `igs_run2.sbatch` (K = 2 round-1 + GAPw), `gapw_only.sbatch`, `gapw_oracle.sbatch`.

- **(b) GAPw** (`fixture_igm/discrete/scripts/dc_gates.py direct` on the round-1 state_eval groups' even/odd maps, K = 2, site mask
  = the sparse-PPCA support): 9 / 9 IgM arms DISCRETE (.03-.27; 3-state arms split in 2: .21-.27), 0 / 9 IgM null and non-arm
  blocks (.47-1.56 or no call).
  - The IgM null's true angle-median split is also not discrete (.47).
  - qtoyflex: closed missed even with TRUE labels (no call: lobes 9 A apart, corridor full of the arm continuum). The arm continuum
    alone split at its true s1 median is called DISCRETE (.13, lobes 30 A apart; reading, not tested: a curved swing leaves the
    chord empty).
  - Only .28 is calibrated (5 % IgM null median - 3 robust sd); the gates are fixed physical or dimensionless constants.
- **Part 4b (p3 knob-free blocks, every block fed):** 5 / 5 arm main blocks DISCRETE (.04-.22), 0 / 3 extras (b5, b7, b8; nearest
  b7 .36), 0 / 13 null blocks (.48-1.65). qtoyflex's closed reader b3 still gets no call (round 1 42 / 58, lobes 3 A).
  The p3 extras are b5, b7, b8 (b6 = a5 main).
- **Part 6 (decision 8 = GAPw for the site layer):** 14 / 14 IgM arm blocks (10 % p2 and p3, 5 % p2) and 0 / 25 false accepts.
  The fission saturation test gets 6 / 14 power with 0 / 25 false accepts.
- **(a)** per-particle held-out tests on sparse-PPCA block coordinates (100 bootstrap replicates, p < .01) fail:
  - 2- vs 1-Gaussian: its Gaussian null sits at -.0025 nats/ptcl (sd .0002), so it accepts anything non-Gaussian; false accepts
    5 / 12.
  - Spikes vs a Gaussian / uniform-segment continuum: misspecified both ways; false accepts 10 / 12, and it misses 10 % a1.
- **(c)** fission map test: not run. hier-ppca measured that it calls adjacent continuum thirds distinct (z 20.5).

**Why:** the coordinator needs a per-site test so the tree accepts discrete site splits and rejects continua (Part 3's seeds cut
the null twin's continuous arms).
**How to apply:** Part 5's fission saturation test fixes both qtoyflex failures without the .28; check it first
([[saturation-fission-site-test-verdict-2026-11-07]]). Use GAPw for IgM-like short-swing sites, with no call = one leaf. Do not trust it on large curved continua
without a node-specific null, nor on minority states adjacent to a continuum. Do not retry per-particle mixture tests at N ~1e5
without a correct continuum model. RNase P's closed corner resembles qtoyflex's closed state next to the arm, so the tree's
own splits must keep carrying it ([[rnasep-closed-state-must-be-retained]]). Related: [[sparse-ppca-pursuit-verdict-2026-11-07]], [[discrete-vs-continuous-site-gap-2026-11-05]],
[[lever3-latent-refs-seeding-verdict-2026-11-07]], [[cross-half-test-needs-band-above-assignment-2026-11-07]],
[[data-fission-merge-test-verdict-2026-11-07]].
