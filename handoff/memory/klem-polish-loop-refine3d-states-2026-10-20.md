---
name: klem-polish-loop-refine3d-states-2026-10-20
description: "BEST so far (synthetic sixc8, 2026-10-20 04:00; split test only after the polish: NMI 0.47/0.48): tree -> final pass -> POLISH = refine3D_states per half seeded with the final states (classification + poses, gold-standard) -> final pass with the re-classified states as candidates: NMI 0.42/0.43 (tree final 0.33/0.35; old true-pose ceiling 0.45/0.42), closed 0.59, arm 0.65"
metadata:
  type: project
---

User's manual-workflow claim confirmed: refinement separates states, via refine3D_states' RE-CLASSIFICATION (prob over
states each iteration), not via better poses. Recipe = oracle_refine.sbatch with STATES = the final pass's children
(per half, all root particles, consensus start, local 20 -> 8 A, maxits 6, ml_reg=no, NP 4), then tree_finalize6.sbatch
CANDA/CANDB = refined states (refstates_nmi2.sbatch extracts them), POSEA/POSEB = refined poses. Generic driver:
scratchpad polish_generic.sh TAG SEEDDIR NAME [STARTDIR].
- sixc8r (6 open + closed 18 at 8 %): final-pass states NMI 0.33/0.35 -> refined states 0.46/0.47 (34 % / 30 % moved,
  closed 18 purity 0.57/0.56) -> final pass on them: all 6 survive, 0.42/0.43, clade 0.47/0.48, 47 separates from 49
  (0.76 vs 0.71), closed 0.59/0.59, closed-arm 0.65. Tree leaves at the polished poses only: 0.36/0.39 -> the gain is
  the classification. Seeded with TRUE labels it ends at 0.52 (attractor at this SNR); poses 3.5 deg either way.
- Pending at write time: 2 % closed (polish_sixc2r), second round (polish2_sixc8r), REAL gs1 polish (9 states,
  polish_rnasep692_gs1, real_polish_chain.sh: proxy tracking + final pass + closed-arm).
- ROUND 2 (seed = round-1 final states, start = round-1 poses): refined states 0.47/0.48, final pass 0.43/0.43, poses
  3.5 deg -> converged after one round. The cross-half EM in the final pass LOWERS the refined states' NMI every time
  (0.46 -> 0.42, 0.47 -> 0.43): it re-assigns at <= 20 A, refine3D_states at 8 A; per-half refined states are already
  independent, so test them without the EM (XEM=0 EXTRA="--softem 0", run final_xem_polish2_noem).
- 2 % CLOSED (polish_sixc2r): refined states 0.40/0.43 (from 0.33/0.33), final pass 0.38/0.40, clade 0.47/0.48, but
  NO closed state: refine3D_states keeps K = the seeds and no seed was closed (18 <= 6 % of any state) -> it sharpens
  existing states, cannot create a minority one. Real-data answer: seed the closed group explicitly (closed_seeds.sbatch:
  gs1 final states + f692e cand label 64 as state 10; run polishcl_rnasep692_gs1, real_closed_polish_chain.sh).
- SPLIT TEST ONLY (final_xem_polish2_noem: XEM=0, EXTRA="--softem 0"): all 6 polished states pass (pair stats
  0.82-0.96 vs threshold 0.125) -> final NMI 0.47/0.48, clade 0.53/0.54, closed 0.56/0.57 = BEST. Recommended order:
  tree (proposals) -> cross-half final pass (K) -> per-half refine3D_states seeded with the final states (1 round is
  enough) -> split test only (no EM). Oracle-seeded refinement attractor 0.52 = the ceiling of this recipe.
- 2 % with an ORACLE closed seed (polishcl_sixc2r: final states + the 1000 true closed per half as state 7): the
  refinement keeps a closed-ish state but it grows to ~4.1k at 18-21 % closed (absorbs 42/46, the nearest open
  conformers); final pass keeps 7 states, that one closed-arm 0.53 -> at 1000/half per-particle evidence cannot hold a
  pure closed state even when seeded (= the real f692e unprotected behaviour). Real data has ~7k closed per half.
- REAL plain polish (polish_rnasep692_gs1, gs1's 9 final states, NP 10, ~55-70 min per half): 52 % change state; poses
  move median 3.7 deg but p90 20.4-20.7 deg (synthetic polish p90 4.5; per-node real refinements p90 5-6.5) -> a large-
  move tail on real data, unexplained; final pass keeps all 9; known-closed proxy stays 5 % of state 2 (as gs1), max
  closed-arm 0.54 -> no closed seed, no closed state (as on synthetic 2 %).
  Tail anatomy (pose_change_hist.py): 7.6 % rotate 135-180 deg (flips); by state: 6 (the gs1 junk-like state, consensus
  corr 0.86, 0 % closed proxy) 57 % > 20 deg (median 43), states 4/8/9 26-34 %, states 1/2 4-5 % -> per-state pose
  instability under refinement flags junk / ambiguous-pose particles.
- FULL 55-conformer landscape (dtree_55_200k_lat, ~98k per half, junk): polish (10 seeds) -> refined states clade NMI
  0.30/0.32 (conformer 0.12/0.13), final pass with EM 0.28/0.28 vs the tree's final 0.25/0.26 -> the gain generalises,
  smaller on a continuum; junk stays isolated (0.90-0.97 in junk-majority states). (polish_generic2.sh MODE arg;
  its refstates readout used the 55_700k truth -> ignore those lines, rerun refstates_nmi2 with MODE=55_200k.)
- REAL CLOSED-SEEDED polish (polishcl_rnasep692_gs1: gs1's 9 states + f692e group 64 as state 10): the seed grows
  2.6k -> 34.5k / 33.4k but holds 74 % / 73 % of all known-closed proxies at 11.2 % / 12.4 % (trees: 5-6 %); all 10
  pass the split test (0 merges); after the EM: state 5 (40k) 77 % of proxies at 9.9 % / 11.2 %, closed-arm 0.69 (best
  whole-particle map; f692e frozen group 0.79). Poses p90 13.5 deg. -> seeded + validated, but diluted ~13x: the arm
  is a small part of the particle; refine3D_states has no focus mask (only klem boards: focal mask) -> next: focused
  (arm-masked) classification per half, seeded the same way.
Related: [[klem-dtree-per-node-refinement-2026-10-20]], [[klem-data-derived-schedule-2026-10-20]], [[klem-board-rewrites-shared-project-trap]].
