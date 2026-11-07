---
name: flexible-target-shortlist-empiar-2026-10-22
description: "lit search 2026-10-22: next real targets for the state tree after RNase P, with EMPIAR IDs, what is deposited and what the source papers left unresolved; OLE RNA 12707 and telomerase 11985 ship particle stacks; CAHRA 2026 benchmark sets 13973/14006/13991"
metadata:
  node_type: memory
  type: reference
  originSessionId: 5f12aee4-62b6-4518-9284-7a4df90b8e47
  modified: 2026-10-23T03:06:52.607Z
---

The user asked for flexible molecules on EMPIAR where heterogeneity methods have not succeeded. Nothing was downloaded
([[no-empiar-downloads-synthetic-igm-pdh]]).
- OLE RNA dimer EMPIAR-12707: 87,716 particles with poses (437 GB). Kretsch Nature 2025: P9.3 flexible, other conserved parts
  unresolved, no heterogeneity analysis tried. ROOL 12708 (31.6k ptcls) and GOLLD 12709 (20.4k) have disordered internal
  linkers. HEARO, the paper's flexible failure case, is NOT deposited.
- Human telomerase: EMPIAR-10732 (43,639 movies, Ghanim Nature 2021) and 11985 (particle stacks 144k-199k). Lobes are tethered
  by RNA. Multibody refinement showed multiple lobe conformations, with the core at 9-10 A; per-lobe maps 3.4/3.8 A by
  subtraction.
- Pks13 EMPIAR-11608 (Kim NSMB 2023): ACP1 in two positions 60 A apart (3-4:1); ACP2, ACP3 and TE never resolved (~20 A
  classes); only focused 3D classification was used.
- IgM 11077 and PDH 10489 are in [[igm-pdh-papers-heterogeneity-baseline]]. Native C. thermophilum OGDHc 11502 and PDHc 10625
  are alternatives to PDH. Dynein-LIS1+ATP 12715 (Nguyen Nat Commun 2025) is a discrete+flexible mixture.
- CAHRA 2026 challenge (Feathers, Zhong et al. bioRxiv 10.64898/2026.09.15.751515), with results out: 13973 is an experimental
  compositional mixture with reference labels; 14006 is a simulated conformational set (1.75 M particles, ground truth); 13991
  is a simulated pose-entanglement set (30k).
- Tetrahymena ribozyme is KEPT; the user challenged dropping it for "2.2 A = rigid", and a global resolution is no rigidity
  test (RNase P's core was < 2.6 A too). Misfolded set 11267 (PNAS 2022): 42,382 micrographs, 10.3 M picks, 2.05 M consensus;
  heterogeneous refinement only gave M1 40.8 / M2 21.0 / M3 8.9 / N 29.3 %; the M core is poorly resolved and P7 flexible.
  In 11844 the 1.18 M deposited particles = exactly the two final maps (473k + 708k), i.e. pre-selected; re-pick the 18,366 raw
  movies instead. Su 2021 (not on EMPIAR): 25 % of picks kept, P9.2/P13 at 4.6 A.
- Dropped: the group II intron RNP 11315/11316 (rigidified by the maturase).
  Integrin avb8 10345 is already a 3DFlex/GaussianEM benchmark, so it is a comparison set, not an unsolved one.
