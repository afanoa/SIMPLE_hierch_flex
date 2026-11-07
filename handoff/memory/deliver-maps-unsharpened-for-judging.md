---
name: deliver-maps-unsharpened-for-judging
description: the user judges state maps at ~10 A without a B-factor (what a class average shows); B -150 at 7 A hid the closed arm the map actually contained; deliver lp10 unsharpened + mild B-50 lp8 versions, not only the sharpened ones
metadata:
  type: feedback
---

2026-10-16 ~01:40: the within-10k docked RNase P map (7.1 A, arm occupancy 0.97) delivered with B -150 at its FSC
looked "closer but a long way to go" to the user; the same map at 10 A unsharpened and at 8 A with B -50 was
"yup those latest two maps are it". Projections of the map at the poses of the sub-class average the user pointed to
(sub 1413) already matched it at 0.85-0.87 in the arm region, better than the deposited closed map (0.77).

**Why:** a global B-factor at the body's resolution amplifies the body's fine structure and suppresses a peripheral
region whose local resolution is lower; the user's reference for "what the state looks like" is a class average,
effectively ~10-12 A and unsharpened.

**How to apply:** for every delivered state map give the lp10 unsharpened version (and a mild B -50 lp8) alongside or
instead of the B -150 one; before concluding a feature is missing from a map, render the map as projections at the
poses of the average that shows the feature (render_subclass.py). Related: [[rnasep-klem-real-data-setup]],
[[rnasep-klem-mask-clips-closed-arm]].

EXCEPTION (user 2026-11-03, RELION leaf refines at ~5 Å): "I want sharpened maps not low pass filtered". For high-resolution
RELION refines deliver PostProcess/<job>/postprocess.mrc (sharpened, FSC-weighted, unmasked), never run_class001.mrc (unfiltered
to Nyquist, looks noisy) and never a home-made low-pass.
UPDATE same day: the global-B postprocess.mrc maps were rejected ("they suck"): noisy solvent, fragmented periphery. Do not deliver
global-B sharpened maps; the candidate is relion --locres filtered maps (job 54931562). Stop any running copy the moment a map type is rejected.
FINAL (same day): after global-B, RELION locres and SIMPLE nu_filt3D versions, the user said "Raw looks the best". For RELION
refines deliver the raw run_class001.mrc only; no filtered/sharpened variants unless asked.
RNase P 2026-11-05: CTF-refined per-leaf refines (Refine3D/leaf_<L>, "ctfref") also rejected next to raw leafd_<L>/run_class001.mrc; deliver and render leafd only. See [[rnasep-relion-polish-pilot-2026-11-03]].
