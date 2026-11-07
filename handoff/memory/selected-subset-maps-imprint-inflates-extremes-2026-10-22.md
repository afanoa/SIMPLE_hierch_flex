---
name: selected-subset-maps-imprint-inflates-extremes-2026-10-22
description: "TRAP: a map of a data-selected particle subset rebuilt from the pixels that selected it overstates its extremeness (swing30 truth: +4-7 A per end; f692st r1 closed-arm own-pixel ~0.2 higher); judge extremes on held-out-pixel maps; f692st closed leaf r1_s5 re-assigned held-out 0.78 (own 0.97-0.99, tree map 1.00)"
metadata:
  node_type: memory
  type: project
  originSessionId: 5e5b4054-6080-49b8-8249-8b67e5b28864
  modified: 2026-10-23T19:06:43.863Z
---

Measured 2026-10-22 with the landscape stage's pixel split ([[landscape-sampling-stage-2026-10-21]]): Fourier pixels split by the parity
of kx+ky on UNMASKED images carry independent noise, so a subset selected with pixel set A can be rebuilt from set B (lsamp2h.py writes
heldout_X / own_X / final_all maps; lsamp_swingpos.py places maps on the true swing direction conf200->conf230).
- s30a r22, the same fold cells: held-out maps 21.1 -> 8.7 A vs their particles' TRUE mean swing 22.8 -> 6.8 A; own-pixel maps 29.4 -> 0.2;
  all-pixel final samples 27.3 -> 2.3 (truth 23.5 -> 6.1). Own-pixel readouts exaggerate each extreme by 4-7 A; held-out ones track the
  truth (slightly compressed at the ends: half-pixel coordinate).
- Real f692st r1 (1-D landscape, 7 samples): closed-arm own-pixel 1.17 -> 0.20, crop all-pixel 1.09 -> 0.28, SIMPLE viewing maps
  1.22 -> 0.39, held-out 0.86-0.87 -> 0.46. f692rf r22 grid closed end: SIMPLE 0.76, crop all-pixel 0.69, held-out 0.55-0.58.
- The tree's own leaves (lsamp2h.py --leaf_heldout 1: re-assign the node's particles to the tree's leaf maps on one pixel set,
  rebuild from the other; agreement with the tree's labels 57 % r1 / 65 % r51): f692st r1_s5 re-assigned 13.3-13.6k, own-pixel
  0.97-0.99, HELD-OUT 0.78; the landscape's 3 most closed samples at the same N (13.2k) 0.77-0.78, its ~2k end sample 0.86-0.87;
  r51_s5 held-out 0.65-0.67 (~5.6k) vs landscape end 0.62-0.63. Honest readouts: landscape ~= tree's closed leaves.
**TRAP (2026-10-23, full692 conformer-43 test): the kx+ky pixel split LEAKS on the real RNase P stacks.** They are the
REGISTERED phase-flipped raw particles that denoise_project writes (22_denoise_project/denoise_raw_particles_part01-10; the
simrnasep692 caches).
- Each particle was rotated and shifted onto its 2D class frame. That interpolation couples neighbouring Fourier pixels, which
  have opposite parity.
- NOT neighbour-averaging denoising: I first claimed that wrongly; the coordinator corrected it from the denoise_project UI text.
- A phase-randomised FAKE arm template gave held-out (pick A, map B) maps that reproduce the fake arm: CC 0.48-0.70, against
  0.55-0.71 for the real conf43 template and about 0 for random picks.
- It kept 0.29-0.33 of its own-pixel amplitude on held-out pixels (conf43 0.40).
- So held-out readouts on these stacks carry about 31 % of the imprint. That includes the r1_s5 0.78 and landscape numbers above,
  and the per-site survivals.
- Always run a fake-template null first, and require about 0 before trusting a held-out map.
- Raw boxes are not a direct fix: the tree's poses are in the REGISTERED frame. Boxes cut from the motion-corrected micrographs at
  the tree's xpos / ypos score -0.03..+0.00 against consensus projections at those poses; the registered images score +0.24.
- The tree's own unregistered stream stacks (/mnt/beegfs/lea/...) are gone, and reproc/13_extract holds different picks.
- No even / odd frame-sum micrographs exist.
- USE A GAPPED FREQUENCY-BAND SPLIT on registered stacks: pick on shells 14-23, rebuild from 2-10, leave 11-13 out.
  - Validated 2026-10-23: fakes Spearman 0.02 between bands, held-out about 0. Scripts: full692_analysis/sw_bandscore.py and
    sw_bandheld2.py.
  - Test it with a band-COHERENT fake (e.g. a template mirrored through the box centre). A phase-randomised fake reads 0 by
    construction on a band split.
  - Use 20+ random draws for the null: the low band's map CC is noisy (+-0.1).
**Why:** selecting particles by their own noisy data puts that noise into their map along the selection direction; both half-maps
share it, so FSC / half-map agreement cannot reveal it.
**How to apply:** judge extremes (closed-arm occupancy, swing, ends of a landscape) on maps from pixels or data the selection did not
use; do not rank a landscape sample against a tree leaf on own-pixel readouts (the leaf r1_s5 "1.00" is own-pixel too - an honest
comparison needs the leaf's assignment repeated on one pixel set). Related: [[split-test-holdout-and-drain-traps-2026-10-21]].
