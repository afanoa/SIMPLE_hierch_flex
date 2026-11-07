---
name: rnasep-f692-landscape-and-family2-defocus-2026-10-16
description: first full klem tree on the 692k apo RNase P (2026-10-16): 47/62/104/180 states at similarity 0.6/0.7/0.8/0.9, 8 of 10 deposited found, 8 embedding islands; family 2 (37 %, "three-pronged" maps) matches 0-4 of all 55 deposited apo conformers AND is the low-defocus half (parent z1 vs defocus Spearman -0.61) -> suspected contrast/defocus split, pscale test pending
metadata:
  type: project
---

Run f692 (chain5, Fortran distributed boards, 10 workers; family stage 14:30-14:33, regions 14:33-17:49, ALL DONE
17:54; ~3.5 h daytime). Landscape over the leaves: 47 states at similarity 0.6 (median 46k ptcls), 62 at 0.7, 104 at
0.8, 180 at 0.9; deposited conformers identified 8/10 at every threshold. Whole-tree embedding: 8 DBSCAN islands,
all particles in islands (figs/fig_global_embedding_f692.png).

All 55 deposited apo maps now docked (routes/toy55/maps160/docked_all, rnasep_emdb_all.csv: 76 entries, 55 apo + 21
holo, EMD-74070..74300). Family 1 refined states name 10 of 55 (mean 0.54); family 2: 2/55 (0.20), its regions 4/55
and 0/55. Family 2 = 256k = 37 %: refined states 9-16 A vs 6-11 A in family 1, user sees "three-pronged" maps.
CONFOUND: parent-board z1 correlates with defocus (Spearman -0.61); per-cell mean defocus 0.85-1.77 um; family 2 mean
defocus 1.24 vs 1.51 um; within the families the leakage is small (|rho| 0.12). The Fortran board has pnoise but NO
per-particle contrast scale (klem.py --pscale exists). Next: pscale in the Fortran board, rerun the parent, check
rho(z1, defocus) and whether an unmatched family survives, before any claim of a new conformer.

Readout bug fixed 18:35: deposited_match.py scored a single state against itself (0.00 everywhere at level 2); now
deviation from the consensus when K == 1. Toy integration tests all negative: global flex prior (15-NN 0.387 -> 0.376),
cell-weighted PPCA = unweighted cell PCA, node-local flex worse than global (kNN 0.19 vs 0.27). Related:
[[klem-f692-fortran-chain]], [[klem-states-need-data-support-not-counts]], [[python-latent-needs-ctf-weighted-residual]].

Update 18:45: pscale (per-particle contrast scale, now in both Fortran board paths, fboard PSCALE=1) does NOT remove
the split: rho(z1, defocus) -0.61 -> -0.43, junk share 17 -> 21 % and junk-defocus rho -0.83, families 247k/279k at
1.31 vs 1.51 um (identification of the parent cells improved: 7/10 at 0.40 vs 5/10 at 0.33). Next hypothesis: CTF SIGN.
The 10 stacks are the phase-flipped denoised stacks (ctf=flip strings) while the board applies the signed CTF
(ctfabs 0, chosen by a correlation calibration at 16-12 A); at 1 um the first CTF zero is at 14 A, inside the parent's
20-10 A bands, so a sign mismatch grows with lower defocus. Test running: parent with ctfabs 1 (fboard CTFABS=1),
with and without pscale (f692ca, f692caps). If it holds, every chain must run CTFABS=1 and the real landscape is void.

