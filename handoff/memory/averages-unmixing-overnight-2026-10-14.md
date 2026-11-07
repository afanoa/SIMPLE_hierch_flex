---
name: averages-unmixing-overnight-2026-10-14
description: "User 2026-10-14 ~01:20: get discrete states from the EXPANDED CLASS AVERAGES with flex working overnight; NO region-focused flex (rejected: 'why are we doing this focused stuff'); deconvolution/unmixing of the averages is the line to develop"
metadata:
  type: feedback
---

User rejected region-focused flex (pcg_mskfile region masks) outright: "Didnt want to pursue that. i simply want a way to
be able to use the expanded class averages with the flex" and "get this to work overnight. no flex with regions".
They liked the averages-derived motion basis ("better than the 3D latent on particles; continue this line").

**Why:** their goal is states built from the cls_expansion sub-class averages through flex, not new support schemes.
**How to apply:** develop deconvolution (unmixing) of the averages in a latent not dominated by member count; fix flex's
treatment of averages with unequal members; validate on toys (20 deg: pooled-latent unmixing 0.984, ceiling 0.985;
flex-on-averages latent fails, R^2 0.20 vs size 0.42), then RNase P. Related: [[discrete-states-fast-toy-loop-2026-10-13]].

**Results ~02:00 (Python prototype on Fortran latents):** flex's OWN latent of the averages never carries the state
(20 deg: R^2 0.20 vs member count 0.42; 13 deg -0.34; de-/re-shrinking by members does not help) -- the averages'
basis is spent on size/sharpness. Sub-class POOLED particle latent + per-parent centring + weighted archetypal
analysis works: 20 deg K=2 corr(w,f) 0.970, particles 0.984 (ceiling 0.985); 13 deg K=2 0.937, 0.951 (0.957);
10 deg needs K=3 (minority = small archetype 4.5 %, 0.892 of ceiling 0.927); 6/4 deg nothing (sub-classes not
enriched). Without parent centring the extremes are view offsets. Extra K splits the majority on a non-state axis
(13 deg K=3: 67/23, both 5 % minority). RNase P 16-per-class K=4 parent-centred: 65 / 34 % view-balanced (TV
0.04/0.08) reproducible (0.85/0.92) + two 0.5 % S3/S5 outliers; maps (rnasep_coupled/k4maps, ~/Downloads/
rnasep_k4_states_20261014) state3-state2 difference r 0.65 with (S2+S4)/2 - S1 at 15 A (S3/S5 0.17). Hard-split maps
are mixtures (A = mean archetype weights); X = A^-1 M deconvolves (condition 3.9).

**Fortran state_placement=unmix (snap reseed27, ~/SIMPLE_master1013 uncommitted; new src/main/flex/states/
simple_flex_pca_unmix.f90 -> batch2 build needs `cmake .` first; infile2=cls_expansion_labels.txt):** toys 20 deg
0.983, 13 deg 0.948 (Python 0.984 / 0.951). E-step seed MUST be the HARD labels (estep_loop8 SEEDHARD=1): 13 deg
0.9931 -> 0.9940 over 5 its, fractions 0.131 (truth .13); the soft proportions (minority share inflated to 19 %)
drift to 0.877. 10 deg needs K=3/4 + a merge (E-step with K>2 keeps 3-4 states; open). RNase P: K=4 is
START-DEPENDENT (furthest-first 65/34 vs k-means++ 8 % S1x1.85 + 62-68 % + 23-28 %, residuals within 4 %); robust
across both expansions and every K = an S1 <-> S2/S4 axis (e16 K=3: 23.5 % S1 x1.57; e100 K=3: 49/49 %, TV 0.025,
S1 x1.29 vs S2/S4 x1.14-1.15; e100 K=4: 22 % S1 x1.52 / 57 % / 19 % S2 x1.44) plus a 1-2 % small-sub-class junk
group (S3/S5). Runs: rnasep_coupled/unmixplace/e{16,100}_K{2..5}; E-step from hard labels: rnasep_coupled/
estep_unmix/{e16_K3,e100_K4} (estep_rnasep.sbatch, score_states nthr 32, 3 its).

**RNase P end-to-end (02:20-02:45, reseed27):** unmix hard labels -> per-particle E-step (estep_rnasep.sbatch:
reconstruct3D objfun=euclid sigma_est=global bootstraps sigma2; score_states true cross-half, nthr 32: ~3 min per
692k pass; EM fractions). e100 K=4: S1-rich 18.8 % (S1 x1.67), middle 62.8 %, S2-rich 9.9 % (S2 x1.29, S1 x0.46),
junk 8.5 % (S3 x2.00, S5 x1.84; grows 1.4 -> 8.5 % = the E-step pulls junk into its own state); e16 K=3: S1 state
20.6 % (S1 x1.63), main 71.2 %, junk 8.2 %. Label changes 48 -> 30 -> 19 % (e100), 32 -> 19 -> 12 % (e16). Maps:
S2-rich minus S1-rich vs (S2+S4)/2 - S1 at 50-15 A r 0.70 (placement) -> 0.77 (after 3 E-steps); e16 main - S1 state
0.56 -> 0.67; S3/S5 differences stay 0.08-0.23. Continuation its 4-5: estep_rnasep_cont.sbatch. Figures/maps in
~/Downloads/rnasep_unmix_states_20261014 and ~/Downloads/rnasep_k4_states_20261014.

**Final (02:55):** 5 E-steps converge (label changes -> 7-12 %) to the same 3 states from both seeds: S1-like ~20 %
(S1 x1.6, S2/S4 x0.55), main ~65-70 % (S2/S4 x1.15-1.18), junk ~10 % (S3 x2, S5 x1.9); the K=4 S2-rich extreme fades
into the main state (14 -> 6 %, S2 x1.40 -> 1.11) = S2/S4 a continuum, S1 distinct (as in the soft-EM memory). Final
maps (e16 K=3 it5): main - S1-like vs (S2+S4)/2 - S1 r 0.72 at 15 A (e100 middle - S1-rich 0.67). The S1-like state's
S3/S5 enrichment rises with iterations (1.08 -> 1.50/1.57): open. Maps/figure ~/Downloads/rnasep_unmix_states_20261014.
Nothing committed or pushed.

**2026-10-14 later:** origin/master now has cls_expansion subcls_avg=kernel (Hans 6d56e6da1; stream p07 uses it). Built
clean at batch2 ~/SIMPLE_kernel (worktree, snap 5e1e6aeee). User: "I will judge the class averages we get out of this
before deciding on building a flex for the class averages" -> do NOT build the flex class-average mode until they
decide; deliver montages of kernel vs model sub-class averages (toys + RNase P classes 486 477 350 187 in
rnasep_coupled/kernel4/{kernel,model}).

**User verdict 2026-10-14 ~10:40:** "kernel path works fine on the rnasep; we cant do anything with the model route
because of weird contrast" -> use subcls_avg=kernel for RNase P expansions. Measured on 4 classes (486 477 350 187,
16/class): kernel pops 457-1164 (model 35-4374), S1 share per sub-class 0.04-0.59 (spread 39-76x binomial; model
3.5-5.8x), sibling even-odd noise corr +0.14..0.21 (not ~0: overlapping kernels share members). Toys at 40 2D
classes: kernel sub-classes NOT state-enriched (view spread dominates coarse classes) -> rerun with 200 2D classes.
Running: RNase P all 112 classes x 100 = ncls 11200 kernel (rnasep_coupled/kernel100); toys fine2D_200 (8/class).
