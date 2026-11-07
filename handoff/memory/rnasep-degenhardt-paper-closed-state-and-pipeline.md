---
name: rnasep-degenhardt-paper-closed-state-and-pipeline
description: "Degenhardt et al. Nature 2026 (RNA catalysis emerges from dynamic structural ensembles, gstRPR RNase P RNA) - closed = TL-TLR docking of P10.1 onto P12 (IDTM folded), ~13 % of apo; their pipeline (21.5 M picks -> ~50 best/mic, rigid-core alignment, 3D classification with focal masks, NU-refine, >= 50k particles per conformer, box 480 @ 0.732 A = 351 A); our box 234 A clips the arm tip's delocalised signal"
metadata:
  node_type: memory
  type: reference
  originSessionId: a718a461-105d-44de-b65c-c91952cccfae
  modified: 2026-10-19T16:43:47.339Z
---

Paper: https://www.nature.com/articles/s41586-026-11140-z (Degenhardt, Stagno, Wang et al., 9 Oct 2026). 76 conformers
(55 apo, 21 holo) of G. stearothermophilus RNase P RNA. CLOSED = the tetraloop-tetraloop-receptor (TL-TLR)
interaction between P10.1 and P12 in the S-domain formed, with P12's interdigitated T-loop motif (IDTM) properly
folded; OPEN = no TL-TLR, IDTM quasi-helical. Closed ~13 % of the apo ensemble, 43 % with rnpA (holo). Flexibility
within every conformer sits in P10.1, P12 and P19 (periphery, ~5 A local resolution vs < 2.6 A core).

Pipeline (EMDB EMD-74209 header + peer-review file): TFS Krios, Falcon IV, Selectris X, 300 kV, 78,777 micrographs,
0.5-2.0 um defocus; 21,475,910 particles picked (~250/mic), iterative cleaning = 2D classification + heterogeneous
refinement with decoy classes, ~50 best particles per micrograph kept; poses aligned on the RIGID CORE, then 3D
classification (cryoSPARC 4.6) with full solvent masks and FOCAL masks on the peripheral regions, iterated with
particle pruning on subsets; 3DVA/3DFlex ONLY on stacks giving incomplete volumes (ED Fig 2f), illustration only, none
of the 76 maps came from them (Supp Video 1: 3DFlex, rigidity 1, box resampled to 128, open state periphery; Supp Video 2:
3DVA at 6 A filter, one component = TL-TLR open->closed; the paper never says why they weren't the main route); per-conformer NU-refinement +
reference-based motion correction; conformer 18 (EMD-74209, closed) = 66,142 particles at 2.9 A, box 480 @ 0.732 A
(351 A); EMCrafter simulations: >= 50k particles for 3 A; reproducibility by re-derivation (CC 0.99) and reduced
datasets; rare extended conformers deliberately enriched (true populations << 6 %).

OUR data: the same collection (FoilHole_..._20241122 EER, Falcon IV) at batch2 /mnt/beegfs/elmlund/rnasep: 31k
micrographs motion-corrected at 1.464 A/px (2048^2, preprocessing/motion_correct), picks in reproc/12_pick (.box),
extraction reproc/13_extract at box 160 (234 A), denoised stacks 22_denoise_project; 692k refined particles. The
closed arm reaches 114 A from the centre (half-box 117 A). The box-delocalisation explanation was TESTED and
REFUTED (job 53902244): the 40k closed particles split into defocus thirds (0.98 / 1.39 / 1.70 um) give tip density
0.60 / 0.65 / 0.67 of the deposited map (base 0.75 / 0.78 / 0.79): the tip does not fall with defocus. A re-extraction
at a larger box is NOT justified by the data (it would still be possible: picks in 12_pick, micrographs at 1.464 A).
The tip is uniformly weaker than the base: mixture at the tip (focal P12/P10.1 boards) or intrinsic flexibility
(pure-toy tip ceiling job 53902293) remain. Related: [[rnasep-klem-real-data-setup]], [[rnasep-klem-mask-clips-closed-arm]].

CORRECTION 2026-10-16 (verified against the Nature text and the 76 EMDB headers): the "13 % closed in apo, 43 % in
holo" are CONFORMER fractions: 0.13 x 55 = 7.15 and 0.43 x 21 = 9.03, so 7 closed apo conformers + 9 closed holo
conformers = 16 closed conformers (user's number). Per-conformer particles (EMDB number_images_used): apo 45,331
(conf 54) to 156,860 (conf 49), median 80,440, conf 18 = 66,142, conf 51 = 82,634; the 55 sum to 4.61 M. All 55 apo
entries cite 78,777 micrographs; the 21 holo entries (EMD-74265..74289, "holoenzyme in 1 mM Mg2+") are a separate
dataset (998,359-particle consensus). OUR DATA = the single session /mnt/projects/CSB-CryoEM-raw/krios/2024/nov/
20241122_111749_4_apo_RNAseP_grid5 (18,777 movies over 2024-11-22/23: 6,280 + 12,495 micrographs; the earlier "31k
micrographs" counted five files per micrograph), i.e. 24 % of their apo micrographs: expected in our 692k = ~16k
particles of conf 18, ~19k of a median conformer, ~90k closed particles (13 %) over 7 closed conformers. The 10k
within-docked set is the purity core of that ~90k, NOT a ceiling. Other sessions on the raw share: apo 2024-12-23
(31,428 movies), 2024-11-22 apo (11,360), apo gri10 (9,644), grid3 (10,056) -> with grid5 ~81k = the paper's apo set;
holo: 2025-08-04 holo-gst-RNaseP_1mM-Mg (32,951 movies, the paper's holo condition), 2025-07-25 holo 2 mM Mg (29,456),
2025-07-18 holo + tRNA (41,571), 2026-07-31 holo + 3x tRNA (69,914). Nothing of the holo is in the 692k.

METHODS PDF (user-supplied 2026-10-16, Downloads/41586_2026_11140_MOESM1_ESM.pdf, 64 pp): holo sample = purified RNA +
recombinant rnpA at a 1:1 MOLAR RATIO, 10 min on ice, then concentrated (so partial occupancy is possible; nothing in
the methods classifies unbound RNA out of the holo data or rnpA-bound particles out of the apo data; the two datasets
were processed separately with the same workflow: iterative 2D/ab initio/heterogeneous-refinement cleaning with
decoys, consensus + LOCAL refinement on the rigid core, 3D classification with the solvent mask, focal masks on the
flexible regions, volumes < ~35k particles lose resolution). Holo: a ~530k-particle subset has NO S-domain density
(ordered C-domain only) = a disordered-arm population; the 50 %-data test gave 8 holo volumes (5 open, 3 closed) and
28 apo (27 open, 1 closed). Apo/holo in one dataset is therefore not something the paper did; it would arise as
partial occupancy in a holo sample. Holo + pre-tRNA sample: 1:2 molar ratio at 2 mM Mg.
