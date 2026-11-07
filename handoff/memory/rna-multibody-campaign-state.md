---
name: rna-multibody-campaign-state
description: "Where the RNA large-motion campaign stands as of 2026-08-13 evening: pilot 36 differentiating, production 37 auto-queued, NU+pc domain masks shipped, stage-3 per-domain pose prototype is the pickup work"
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T19:17:35.037Z
---

State at 2026-08-13 ~15:30 (user disconnecting for hours; box works autonomously):

**Amplitude ladder measured** (rel diff energy end-vs-end, same metric): flex sextiles 0.47 ->
8 curve frames 0.60 -> 12-frame 4%-tails 0.79 -> refine3D pilot iter9 0.65 with conc 0.33 (flex
conc ~0.23-0.25 = fade; pilot's higher conc = structured displacement). Latent is attenuation-
limited, not saturated; image-domain assignment is the amplitude-faithful path. kd ladder math:
IgG swing z=2*pi*d/lambda~4 -> 3 eigenvols suffice; RNA ~40A at lp=12 -> z~21 -> 1-axis frames can
only breathe (user confirmed the 3-eigenvol IgG rotation).

**Running on fwl-c142339:** `36_refine3D_pilot` = seeded 3-state refine3D_multi
(multivol_mode=input_oris_refine, seeds = 34's endpoint frames + consensus, update_frac=0.25,
staged workflow, state overlap fell 0.77->0.48 = states differentiating, populations 123k/101k/...).
**Auto-queued:** `queue_37.sh` (watcher PID 1191450) waits for 36's master then launches
`37_refine3D_seeded` = FULL 932k 3-state refine3D_multi seeded from 36's finals (falls back to
latest recvol_state0N_iter*_lp if rec_final missing). Logs: 36/nohup.out, 37_queue.log,
37_refine3D_seeded/nohup.out. Est 3-5 h.

**On reconnect, check:** 36's rec_final_state01-vs-03 visually (displacement vs fade -- THE
verdict); 37's progress; then per-domain angle work. Masks for the user sit in
`reextract/nu_domains/` (domain_mask_001 = swept volume from pc lobes, 57A equiv reaching 77A;
domain_mask_002 = nu+pc joint, 36A at r=44-66; core 73% of NU-automask envelope).
`scripts/flex_nu_domains.py` is the generator (NU automask = envelope; Otsu over quantized locres
tiers; pc lobes ranked WITHIN the outside-envelope shell; closing bridges anchors to lobes).
NU-only partition is 29 fragments / 0 alignable -- discrete all-NU multibody is dead; the agreed
direction is the NU-weighted deformation FIELD (sigma(x) from locres tiers, core pinned,
register-then-PCA composition; user approved "that seems interesting").

**Pickup work (local):** stage-3 prototype = per-domain local pose search on subtracted residuals
(core-projection subtraction, small angular bank around consensus pose), validated by correlating
per-particle domain angle against z2 from 32's cache; then the field-model design doc. Curve
placement (project_onto_target_polyline + scripts/flex_curve_targets.py + SIMPLE_COV_TARGETFILE +
SIMPLE_COV_ZSHRINK=0) is wired, IgG-smoke-tested, and validated on RNA (run 35).
See [[flex-pca-must-be-knob-free]], [[flex-pca-solvent-mask-dominates-rna-basis]],
[[flex-pca-gmm-undoes-equal-mass-placement]].
