---
name: flex3d-soft-em-toy-results-2026-10-05
description: "Mac numpy/numba toy of iterated soft EM over flex 3D states: cross-half Wiener refs DRAIN small states even from the truth; consensus+Wiener-deviation refs + fixed fractions + consensus+deviation M-step work; gain over one-shot M-est needs good poses"
metadata:
  node_type: memory
  type: project
  originSessionId: 64d2af37-1bd8-4ec4-bcaf-cec5cd36a1e4
  modified: 2026-10-05T03:12:53.731Z
---

Prototype (2026-10-04/05, overnight at the user's request): worktree `.claude/worktrees/suspicious-sammet-9e8bee/scripts/flex3d_softem/`
(softem_core.py numba slice/score/backproject kernels, softem_sim.py toy2 states -> box 64 @3.66 A, CTF, SNR 0.05, consensus
poses from a local search against the true mixture; softem_run.py first driver; softem_em2.py second driver with arm spec
START:REFS:OUT:PI:DAMP:POSE). Data/logs in the session scratchpad (t1 = toyK states 1/7/21/27, t2 = close states 14/15/17/13,
pops .4/.3/.2/.1, N=10k); flex stand-in q0 = calibrated GMM posteriors, acc 0.65. Runs ~35 s per E-step on the Mac.

Findings:
- Cross-half scoring against Wiener-filtered state half maps DRAINS small states into the dominant one (T2 pi13 .104 -> .014
  in 5 its); it happens in ONE E-step even with refs built from the TRUE labels (.104 -> .06). Predictive-variance likelihood
  (unfiltered refs + ref noise + logdet) drains faster. Self-scored (own half) looks stable but is self-inclusion overfitting
  (true-label refs give acc .995 vs .84 honest) and locks in its start labels.
- Fix that works: E-step refs = consensus half map + state deviation, each Wiener-filtered by its own half-map FSC, scored
  cross-half; state fractions FIXED to the flex fractions (pi-EM is unstable both ways); M-step output = consensus +
  deviation filtered by 2FSC_d/(1+FSC_d). The output model alone lifts one-shot M-est a lot (T2 dev corr .25/.52/.41/.39 ->
  .43/.69/.59/.62). At TRUE poses the loop reaches the true-label ceiling (acc .87, dev .81/.93/.90/.90 vs B .70/.91/.84/.85);
  at biased consensus poses the extra gain is small (centre state only): poses are the bottleneck.
- Per-state pose loop against noisy per-state Wiener refs DEGRADES maps (fits ref noise); exact refs would give 2.9 deg vs 4.6.
- Whole-map FSC vs truth is consensus-dominated; judge states by deviation corr / deviation FSC.
Related: [[flex3d-soft-em-proposal]], [[score-states-likelihood-posteriors-toyk]], [[spliceosome-mest-refine3dauto-2026-10-02]].

**Real code (2026-10-05, user: "branch off and work on the real code"):** local branch `flex-softem` (session worktree
suspicious-sammet-9e8bee; base 72a9bbb65 = 7e5d05bb6 + flexref working tree) adds `simple_softem_weights.f90` (E-step after
each refine3D batch search: cross-half via pftc%set_eo flip, gen_raw_euclid_vals x wsqsum, logsumexp over all in-plane x
assigned direction + 8 nn, chi2 c per batch, fractions fixed to the loaded flex weights, overwrites build%flexw) and
`SIMPLE_SOFTEM_REFS=sx` (consensus + deviation refs in read_mask_filter_reproject_refvols, no eo blend). Env switches only.
batch2 tree ~/SIMPLE_softem (rsync copy, 3-min build). TRAP: any lpset mode (lp on cline, filt_mode uniform/fsc/
nonuniform_lpset, refine3D_states) merges even/odd refs -> no cross-half; and projects with qsys_name=slurm spawn worker
jobs even at nparts=1 (set qsys_name=local for in-job runs). toySE (/mnt/beegfs/elmlund/afan/softem; states 14,15,17,13,
real flex start that missed state 13): loop acc .785 -> .900; FSC 10-6A frozen .735/.926/.798, fsx (refs only)
.761/.930/.816, sx (loop+refs) .786/.936/.823; pose err 3.20/3.10/3.04 deg; pair-difference FSC 20-12 frozen .63/.69/.65,
sx .71/.76/.72. Plain loop with standard refs = frozen. Missed state is drained, not recovered.
**Equal-iteration control (minits=6):** frozen6 FSC10-6 .734/.924/.797, fsx6 (frozen weights + sx refs) .776/.931/.824,
sx .786/.936/.823 => the REFERENCES are the lever; recomputing weights adds ~.01-.02 + classification (.785 -> .900).
**10180 real data (existing 2026-10-02 M-est runs):** refine3D_states' rec_final_*_even/odd are identical copies (use
recvol_stateNN_even/odd_unfil). Deviation-from-consensus even/odd FSC: soft 40-20A .47-.57, 20-12A .06-.14, >12A ~0
(hard ~.15-.19 at 20-12, state 2 .38): the 5 states differ reproducibly only at >= 20 A, so per-state maps cannot gain
specific high-res signal; consensus+deviation maps reach the consensus 7.45 A. Script: /mnt/beegfs/elmlund/afan/softem/shr_real.py.
**In the user's M-est recipe (refine3D_states flex=no flex_weights=yes lpstart=8 lpstop=5, toySE):** r3s FSC10-6 .713/.915/.802,
pair20-12 .606/.679/.638, pose 3.03 deg; r3sx (+SIMPLE_SOFTEM_REFS=sx, raw-halves deviation FSC) .791/.933/.857, .730/.802/.763,
2.62 deg = the biggest gain. 10180 r3s vs r3sx (same flex refit r10180/F, nparts=5) submitted 2026-10-05 06:13 (53668101/2).
**hsx (hard multi-state + sx refs, toySE):** FSC10-6 .787/.936/.824 vs hard .740/.929/.796, pair20-12 .712/.766/.720 vs .644/.693/.661, pose 3.03 vs 3.19: sx refs are a general multi-state reference model.
**Real-data A/B (2026-10-05, user recipe refine3D_states flex_weights=yes, +/- SIMPLE_SOFTEM_REFS=sx):** 10180 per-state
res0.143 8.91/8.11/8.49/8.11/8.63 -> 8.49/8.00/8.36/8.00/8.49, consensus 7.45 -> 7.25; RNase P S4 4.04/4.04/4.11 ->
4.04/4.04/4.04 (0.5: 5.21 -> 5.09 for state 2), consensus 3.90 -> 3.84 (as-is arm reproduces 2026-09-29 exactly).
Even/odd deviation-from-consensus FSC DROPS with sx on real data AND on the toy, where truth-based state differences
improve: that GT-free metric cannot judge this change (it rewards shared reference bias). Missed state: flex K=6 on TRUE
poses isolates it (pure), on refined poses the deconvolved-mixture macro-cluster fold collapses 6 -> 3.
**Master port (2026-10-05):** local branch `state-refs-shrink` (from origin/master 41ec883e4; commits 3b5469c33 code + note),
param `state_refs=independent|consensus`; batch2 ~/SIMPLE_staterefs; fast gate 14/14 pass (31-32 s; untouched master 30.2 s,
budget overrun is node variance). Master toy (hard multi-state): FSC10-6 .746/.927/.796 -> .795/.937/.826, pair20-12
.652/.700/.660 -> .727/.776/.727. Missed state: image-domain over-split (8 slots) best 13 purity .59 (consensus poses) /
.76 (2 deg) / .82 (true); flex K=6 true poses .99. Traps: oriops on an existing project drops its stk segment; set labels
with `selection oritype=ptcl3D infile=<state per line>` after import_particles (which forces state=1).
