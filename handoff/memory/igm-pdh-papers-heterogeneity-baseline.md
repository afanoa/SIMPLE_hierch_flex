---
name: igm-pdh-papers-heterogeneity-baseline
description: "what the source papers of EMPIAR-11077 (IgM, Chen 2022 Nat Commun) and EMPIAR-10489 (N. crassa PDH, Forsberg 2020 Nat Commun) achieved on heterogeneity, and the real motion/occupancy ranges vs our synthetic arms/tethered sets"
metadata:
  node_type: memory
  type: reference
  originSessionId: 778a2ada-672f-41f9-9779-5d5602884567
  modified: 2026-10-19T18:55:30.253Z
---

Read 2026-10-19 (user asked for a workflow/heterogeneity summary).

IgM, Chen et al. Nat Commun 13:6314 (2022), EMPIAR-11077: crYOLO 1.82 M picks, cryoSPARC 3.2, 301k particles, C1
NU-refine 4.32 A. 3DVA (3 comps, bin 2) clustered into 5 refined states: only the subunit-1 F(ab')2 (held by J-chain
hairpin on Cmu2-1) is resolved, pivoting -20 to -45 deg in-plane. Fabs of subunits 2-5 never resolved: smeared,
range estimated from the envelope at +-50-60 deg in-plane and +-70-80 deg out-of-plane. No cryoDRGN/multibody.

PDH, Forsberg et al. Nat Commun 11:4667 (2020), EMPIAR-10489: 4,887 mics, Arctica/Falcon II 1.25 A, RELION 3.0.7.
I map 69,538 ptcls 4.1 A; T-symmetric classification of the interior protein X gives X 18,558 (4.4 A) / Y 21,129
(4.3 A) / N ~43 %. Periphery (tethered E1/E3) never attempted: "position probability distribution", averaged out.

Only success on PDH periphery anywhere: Wang et al. Nat Commun 16:917 (2025), cryo-ET of porcine mitochondria,
symmetry-expanded 60 sites ~20 nm out + RELION classification: E1 1-46/core (avg 21, bimodal), E3 0-12 (avg 4),
E1 radial SD ~24 A. Our synthetic sets are milder: arms +-40 in / +-30 out of plane, tethered 12 A jitter at 50 %
occupancy ([[klem-synthetic-benchmarks-arms-tethered]], [[no-empiar-downloads-synthetic-igm-pdh]]).
