---
name: pptx-qa-quicklook-and-pymol-level-traps
description: "TRAPS for slide work on the Mac: python-pptx speaker notes make macOS Quick Look hang; PyMOL normalize_ccp4_maps (default on) breaks raw-unit contour levels"
metadata:
  node_type: memory
  type: reference
  originSessionId: dc17548d-1b94-4fe8-ad02-93f812b3f2a7
  modified: 2026-11-04T22:33:44.009Z
---

Slide QA without GUI apps (2026-11-04): no soffice anywhere; render each pptx slide with `qlmanage -t -s 1600` on
single-slide copies (python-pptx venv in scratch, drop the other sldIds). Wrap qlmanage in `perl -e 'alarm 60; exec @ARGV'`.
- python-pptx `notes_slide` adds a default notes master that makes Quick Look (Finder preview) HANG on the whole deck;
  ship decks without python-pptx notes. The skill's validate.py needs Python >= 3.10 (Mac has 3.9): run it on batch2
  with module cryodrgn/3.5.3 in a small sbatch.
- PyMOL loads CCP4 maps normalised (normalize_ccp4_maps 1). Levels computed in raw units (mrcio level_for_volume) need
  `cmd.set('normalize_ccp4_maps', 0)`; the OLE deck's cryoDRGN volumes were drawn with raw levels on normalised maps, so too
  low a contour (noisy, oversized). OLE analysis2.json levels are in normalised units.
Related: [[never-launch-gui-apps-on-the-mac]].
