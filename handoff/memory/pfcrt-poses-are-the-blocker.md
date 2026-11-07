---
name: pfcrt-poses-are-the-blocker
description: "2026-09-11 PfCRT verdict: (1) the 'helices move between class averages' premise is VIEW GEOMETRY (5-11 deg tilts explain 92-98% of the reproducible class-pair difference; 2-3 deg pairs differ only by noise); (2) a third of the particles have poses SIMPLE and the deposited RELION refinement (EMPIAR-10330) disagree on by >30 deg, but removing them reproduces run-56 pc1 at cos 0.92 -- poses are noise, not the source of the eigenvolume; (3) 2D classes lump views. Recipes + traps"
metadata: 
  node_type: memory
  type: project
  originSessionId: e5141833-c145-408a-ba48-7c3a3ffe0fbd
  modified: 2026-09-11T18:59:13.171Z
---

**Finding (2026-09-11 afternoon):** on PfCRT (16.9k ptcls, 9_refine3D_auto at 3.93 A) per-particle 3D poses are unreliable
for ~1/3 of the set (this was FIRST thought to explain the non-localised eigenvolumes; the ablation below says it does not):
- deposited star `EMPIAR-10330/data/03_Final_Refined_Particle_Stack/FinalRefinement-OriginalParticles-PfCRT.star`
  (same 16905-particle stack, same defocus order; copied to `/mnt/beegfs/elmlund/pfcrt/empiar_meta/`) vs SIMPLE ptcl3D:
  after a global rotation fit, 53% agree within 10 deg, 67% within 30 deg, p90 = 126 deg. The deposited refinement was
  run on cryoSPARC micelle-SUBTRACTED particles (`rlnImageOriginalName = J407/imported/subtracted.mrcs`).
- 2D class averages rebuilt (make_cavgs with the abinitio2D params) from the "agree <10 deg" subset (8936) reproduce the
  class averages sharply; the "disagree >30 deg" subset (5610) gives a weak, blurred helix bundle. Selection files in
  `~/pcg_tools/cls/sel_{agree,dis,Bagree,Bdis}.txt`, projects `/mnt/beegfs/elmlund/pfcrt/cls_ablation/R{agree,dis}/`.
- SIMPLE's refinement is "converged" on those particles anyway: corr 0.38 for everyone, <0.2 deg movement/iteration --
  the objective is micelle+Fab dominated and cannot discriminate the transporter pose.
- run-56 latents carry no 2D-class information for the inside components (between-class F 1-3) but the complement
  (micelle) block does (F 14-40): the 2D classes are micelle/view driven.

**Trap:** the 2D classes LUMP views. Only ~25-50% of a sharp top-view class's members lie within 15 deg of the class axis
in EITHER refinement, and the off-axis+agree members' 2D average is a blur, so "the helices move between class averages"
can be a view difference (class axes 2/3/8/21/32/43 span 3-16 deg; a 10 deg tilt moves helix tips ~5 A). The
"pose-off vs class axis" split (selA/selB) conflates true off-axis views with mis-posed particles -- use the
two-refinement agreement instead.

**Ablation RESULT (runs 72 = agree 8936, 71 = disagree 5610):** 72 pc1 = 56 pc1 at |cos| 0.92 with the same localisation (top-10% 0.77, patch radius 21 A, centroid 9 A off centre), het-res 23.7 A (fewer particles); 71 gives a different, worse basis (het-res 25.6 A, cos 0.66-0.77 to 56/72). So the pose-ambiguous third is NOISE, not the source of the leading component; the well-posed majority defines pc1.

**Premise test (pairview.py, synthetic reprojections at the members' own 3D poses, class pairs aligned in 2D incl. mirror):** for same-axis class pairs 5-11 deg apart (2-21, 2-8, 8-32, 3-32) the observed class-average difference is reproducible (even/odd 0.68-0.89) and correlates with the VIEW-ONLY synthetic difference at 0.80-0.91, explaining 92-98% of its reproducible energy; for pairs 2-3 deg apart (2-3, 21-32) the observed difference is not reproducible (0.06, 0.42). The helix 'motion' between class averages is projection geometry, not conformational change. Do not use class-average helix differences as evidence that flex is missing a mode.

**Ablation launched:** run-56 config (protein support, COMPLEMENT=4, EM_DEFLATE=4, DEFLATE_BG=1, box 100, budget 0) on
Ragree vs Rdis: scripts `flex_pca_pose_{agree,dis}_script`, outputs `FLEX_PCA_pose_{agree,dis}_OUTPUT`, jobs
53189789/53189790. Readout: eigenvolume localisation (top-10% energy fraction), init-deflated half match, class-F of z.

**Tools:** `~/pcg_tools/cls/` (p2d.txt/p3d.txt = print_project_field dumps, clsview*.py, compare_*.py, vR.npy,
res_simple_vs_relion.txt). make_cavgs needs `nparts`; reproject needs an EVEN `nspace` + `pgrp` even with oritab and
does NOT apply shifts; `selection infile=` needs one 0/1 line per ptcl2D entry. Montages in ~/Downloads/pfcrt_cavgs.

**How to apply:** before any flex battery on a new dataset, check pose reliability against an independent refinement or
a 2D-class consistency test; on PfCRT, fit flex on the agree subset or re-refine with the micelle subtracted / masked
out (what the depositors did). See [[flex-pca-complement-block-and-step0]], [[flex-pca-validate-inputs-before-batteries]].
