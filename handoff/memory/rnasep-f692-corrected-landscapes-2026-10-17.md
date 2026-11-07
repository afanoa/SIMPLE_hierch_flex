---
name: rnasep-f692-corrected-landscapes-2026-10-17
description: corrected 692k landscapes (abs CTF, pscale, support criterion): f692c with family refinement 60/83/110/184 states, 10/10 docked conformers; f692d without 64/86/138/223, 9/10; 8/12 and 10/12 regions (6 boards failed); dev loop on 150k real data 50 min per family
metadata:
  type: project
---

2026-10-17 results (runs /mnt/beegfs/elmlund/afan/rnasep_klem/f692c, f692d; chain8 with node_states2 --support):
- f692c (refine3D_states per family, 20->8 A maxits 6): states at sim 0.6/0.7/0.8/0.9 = 60/83/110/184, median 48k/34k/26k/13k
  particles per state; all 10 docked conformers (1,10,18,27,28,37,46,48,49,51) identified at every threshold; embedding
  fig_global_embedding.png has no defocus structure.
- f692d (no refinement): 64/86/138/223; 9/10 (conformer 46 only at 0.9) and more states per conformer. Verdict: refinement
  keeps the count and cleans the identification (same as the toy A/B).
- Caveat: 6 of 24 region boards failed (worker-script cwd collision, then the 1 h first-reply timeout under the core cap);
  landscapes are over 8/12 and 10/12 regions; the missing regions were resubmitted on the new build at 12:10 for a recount.
- Dev loop (f692dev, 150k subset, chain9 DEV=1 SERVER=1): family 1 done 50 min after launch (child 7, refine 6, feedback 8,
  3 served regions + level 2 27); family 2 refinement took 42 min (10 states, 20 parts, 32 threads) - the variance to study.
- Dev-run refinement (fam2, 75k, 10 states, 20 parts): 42 min, of which ~24 min master-side volume assembly = I/O over parts x states x 2 halves partial volumes (32 threads did not help); fam1's refinement aborted silently (`Input volume recvol_state04.mrc does not exist`, exit 0, no iteration) -> family script now requires an ITERATION line; next lever = supported K (SUPPORTK=1) + 10 parts.
- all-55: 30 of 55 deposited apo conformers are the best match of >= 1 leaf state map, both f692c and f692d (readout55/identification55.txt).
- f692dev (150k, 6 regions, 2 h with queue waits): 38/67/105/187 states, 8/10 at 0.6-0.7, 10/10 from 0.8.
- STRICT recovery (recover55.py: a conformer counts when a state matches it at deviation corr >= 0.4 with a 0.05 margin over every other conformer): f692c leaf states 8 of 55 (3 7 26 28 31 40 47 49); consolidated at 0.7 only 3 (7 26 47) - whole-map merging loses conformers. This is the progress metric toward the paper's 55 (+21 holo).
- all-55 readout = readout55.sbatch (every leaf state map vs maps160/docked_all), results in <run>/readout55/matches.txt.

**Why:** these supersede the VOID signed-CTF landscape ([[rnasep-f692-landscape-and-family2-defocus-2026-10-16]]).
**How to apply:** quote f692c as the current real-data landscape; report region completeness alongside any count.

**f692e (full recipe, full scale, 2026-10-17 18:37-21:13 = 2.6 h wall):** 49 family states (25+24), region boards 12/12
done, landscape 75/109/156/253 cell-states at similarity 0.6-0.9, docked conformers 8/9/10/10. Zoom (zoom_leaves.sbatch
FAM=1/2) and the data-derived merges (pooled_dd.sbatch) launched on it; see [[klem-region-local-pooling-2026-10-17]].

**f692e final sets (2026-10-18 00:xx):** leaves with ALL particles (UNCAP=1, confidence-weighted ownership): 63
states 3-39k ptcls, FSC median 9.4 A (best 5.3, 10 at <= 7 A), 10/10 docked, strict 7/55, closed max 0.78; 19 junk
suspects (no deposited match < 0.25 or FSC > 10.5 A). Zoom with all particles: 212 sub-states <= 4.9k, strict 11/55
(33/55 identified), pure closed = state_01 (4.9k, occupancy 0.94) + 0.87/0.83/0.83. Delivered to
~/Downloads/rnasep_klem_20261015/deposited_conformers/{f692e_leaf_states_allptcls, f692e_closed_substates}.
Direction-bin 2D route experiment on the 4 synthetic sets: synth_bins_route.sbatch + avgs_states3d.py (K-state
cross-half EM over sub-class averages with known directions), jobs 53921610-3.
