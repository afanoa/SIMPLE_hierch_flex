---
name: igm-retain-side-views
description: "STANDING (user 2026-11-02): IgM picking/2D selection must retain side views (thin ~310 A rods), never drop them as filaments"
metadata:
  node_type: memory
  type: feedback
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-02T22:29:25.796Z
---

IgM pentamer (EMPIAR-11077) is a flat star; edge-on views are thin ~310 A rods that look like filaments or ice edges and get discarded by default.

**Why:** user said "we need to retain side views" while building the crYOLO training set from template picks + 2D cleanup.

**How to apply:** in every IgM class selection, list side-view/tilted classes separately and keep them (check rod-like classes against edge-on EMD-13921 reprojections in ds_batch/4_reproject); include side-view averages (or edge-on reprojections) in pickrefs; report the side-view share of crYOLO training boxes and check rods are picked on pilot overlays. Related: [[igm-choose-vs-em-replays-2026-10-28]].
