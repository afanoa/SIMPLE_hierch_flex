---
name: hear-the-target-sampling-not-the-import-number
description: "when the user repeats a pixel-size instruction, find the target sampling they want downstream before arguing; \"import at 0.7118, bin 2x to 1.42\" meant work at 1.42 A"
metadata:
  node_type: memory
  type: feedback
  originSessionId: a9a4000f-ee34-4500-a4e0-aed13b26e342
  modified: 2026-10-19T15:03:10.609Z
---

2026-10-19, EMPIAR-10509 ribosome: the user said "import at .7118" three times. I kept explaining that the
super-res movie files are 0.3559 A (true, and data-verified) and left preprocessing at 0.7118 A micrographs. What
they wanted was "everything downstream binned 2x to 1.42". User: "How many times do I have to tell you?"

**Why:** the user thinks in physical-pixel terms (0.7118 = the K3 physical pixel). Their repeated number was a means;
the goal was the processing sampling. Arguing the means while missing the goal cost an hour of compute and their patience.

**How to apply:** when a pixel-size or sampling instruction looks wrong, first ask yourself (or ask them in one line)
"what sampling should the micrographs and particles end up at?". Implement that target with SIMPLE's true
`smpd` plus `smpd_downscale`, and state the mapping in one sentence (super-res file smpd = half the physical pixel).
Don't relitigate the import value. Related: [[empiar-10509-ribosome-inputs]], [[act-on-the-literal-ask]].
