---
name: rnasep-flex-seeds-vs-refine3d-finals
description: "2026-09-19/20 measured: flex_pca's RNase P seed MAPS correspond to the refine3D_states finals at 2.5-3.3x lower contrast, but the PARTICLE partitions do not overlap (adjusted Rand 0.025 sphere seeds, 0.004 envelope G10): refine3D_states built its states itself from five enriched seed volumes; the envelope basis is the most reproducible (10/10 >= 0.85) yet gives the faintest k-means states; both flex partitions carry a defocus gradient"
metadata:
  type: project
---

Runs: `/mnt/beegfs/elmlund/afan/rnasep_states/25_refine3D_states` (flex init on 932,496 ptcls: gridding,
sphere mskdiam 200, 4 probe iters, k-means 5 states) and `rnasep_flex_state{1,2,4}/1_flex_pca` (per-state,
pcg + automask 7.2-7.7% of box, 4 iters, 4 k-means states). Scripts `~/rnasep_check.py`,
`~/rnasep_check2.py` on batch2 (numpy 1.14 on the login node is fine for box-160 maps). Full write-up:
`doc/implementation_notes/flex_rna_symmetry_multibody_investigation_2026_09_19.md`.

**Numbers (10 A low-pass, 200 A sphere):** pairwise rel. difference energy: refine3D finals 0.114-0.404
(top-1% conc 0.56-0.75); flex seeds of the same run 0.064-0.190 (conc 0.57-0.67); per-state flex k-means
states 0.011-0.103. Deviation energy from the 5-state mean: seeds 0.029-0.069, finals 0.081-0.212. Final
state j descends from seed j. `STATE OVERLAP` 0.332 after refine3D_states' first iteration -> 0.817 at the
end; populations moved 281k/272k/140k/129k/111k -> 220k/284k/90k/238k/101k.

**Why:** the latent is a continuum (GMM AUTO: 4 macro-clusters each "continuum"; deconv signal fraction
0.75-1.02 per axis, so per-particle noise is NOT the issue, unlike IgG). k-means cells on a continuum
average over conformation ranges = "fade". refine3D_states sharpens references and boundaries by
classification EM at full band with pose refinement; flex_pca embeds once and averages, by design.
Basis: 10/10 matched at it=1 (lead |cos| 0.997), cross-half S(t) FALLS with iterations (0.93-0.97 -> 0.68-0.81
at it=4) while -2logL/N still drops -1.1/iter -> more probe iterations would hurt the basis (cnga1/TRPM4
pattern). Merged spectra flat and non-monotone (342-567). PCs: signedness |s|<=0.04 (transport), conc
0.62-0.79, 0.1-0.5% energy outside support. Support: at most 2.0% of any final's above-3sigma density lies
outside another state's automask (max radius 98-106 A) -> the missing solvent prior is not the RNA limiter.

**How to apply:** treat flex-only state maps as seeds; the flex init -> refine3D_states division of labour
is right for discrete states. Untested flex-internal levers: on-axis equal-mass states along the reproducible
PCs with the MASKED basis (never run on RNA; PfCRT precedent 0.98 end to end), an embedding-only iteration
ladder scored on seed contrast (not eigenvolumes). Traps: `refine3D_states flex=yes` forces gridding + sphere
for the init (85% nuisance in the it-1 basis vs 58-61% with the envelope) and copies `pgrp`; `min_state_frac`
default 0.1 makes `nstates=10` impossible ("floors exceed the particle count", `rnasep_states10/...failed_minfrac`;
rerun at 0.05 gave 10 states of 72k-153k). See [[flex-symmetry-two-regimes-not-expansion]],
[[multibody-flex-design-2026-09-19]], [[rnasep-flex-per-state-plan]].


**CORRECTION 2026-09-20 (label level).** `label_confusion.py` (scripts/flex_pca_verify/region_toy) against the `state`
field of 25_refine3D_states: seed labels vs finals ARI 0.025 (only seed 2 -> final 2 at 45% and seed 4 -> final 4 at 39%
against 30/26% chance), G10 (envelope, pcg, 10/10 pairs >= 0.85) vs finals ARI 0.004. The seed->final MAP correspondence
came from the seed volumes, not the labels. G10's k-means states are the faintest of all (pairwise 0.02-0.05, best-match
to finals 0.38 vs 0.63 for the sphere seeds); mean defocus spans 1.21-1.39 um across seed states (nuisance in the
partition). Never score RNA flex runs by "match to the refine3D_states finals"; use the even/odd reproducibility of
on-axis deliveries (`rna_axis_score.py`). Runs: `/mnt/beegfs/elmlund/afan/rnasep_regions/{G10,G30,R_G10cov_*,C_G10cov_B3}`.
