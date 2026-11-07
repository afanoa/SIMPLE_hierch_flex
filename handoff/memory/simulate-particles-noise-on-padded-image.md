---
name: simulate-particles-noise-on-padded-image
description: "TRAP: simulate_particles adds its snr noise to the 2x PADDED image (var over boxpd^2) then clips, so box SNR ~ 4 x snr; own-noise toys (qtoy) must use the padded variance; snr >= 5 adds no noise"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-25T21:33:25.130Z
---

2026-10-25, qtoy (routes/toy55/run_qtoy_sim.sbatch, qtoy_merge.py).
- simimg adds noise only when snr < 5. add_gauran uses noise variance = var(image) / snr over the WHOLE image it is given, and
  simulate_particles gives it build%img_pad (boxpd = 2 x xdimpd) before clipping to the box.
- So the toys' calibrated "snr 0.0012" (sites toy, RNase P absolute noise) is about 0.005 in the box: the padded variance is ~1/4
  of the box variance.
- qtoy simulates noise-free (snr 10) and adds its own fixed-sigma noise. Using the box variance made the noise 4x too strong:
  refine3D from 10-deg perturbed poses diverged (median pose error 57 deg vs the sites toy's 3.5 deg). Fix: variance of the clean
  image zero-padded to 2x the box (vpad in qtoy_merge.py).
**Why:** a toy that adds its own noise must match simulate_particles' effective SNR, or refinement fails and nothing downstream is
comparable.
**How to apply:** for own-noise simulations compute sigma^2 = mean padded-box clean variance / snr; check the refine3D pose error
(pose_error.py) against the sites toy's ~3.5 deg before running trees. Related: [[hub-toys-and-pose-bias]].
