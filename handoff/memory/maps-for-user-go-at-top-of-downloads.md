---
name: maps-for-user-go-at-top-of-downloads
description: "maps the user should look at go directly in ~/Downloads (top level), not nested in subfolders; halves/README can stay in a subfolder"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-23T19:17:29.038Z
---

When the user asks to see a map ("download the map for it and I will have a look"), put the full map(s) at the TOP level of
~/Downloads. 2026-10-23 I put conformer 43 and then the RNase P junk-group maps in subfolders
(deposited_conformers/, rnasep_full692_final/junk_group_maps/); the user said twice "please put at top of downloads" and
"put these at the top of downloads im not looking for them".

**Why:** the user opens maps straight from Downloads and won't dig through folders.

**How to apply:** copy the full maps to ~/Downloads/<descriptive_name>.mrc (cp -n, never overwrite). Keep halves, READMEs and work
files in a subfolder if needed ([[half-maps-always-in-halves-subfolder]]), deliver maps unsharpened ([[deliver-maps-unsharpened-for-judging]]),
and name the top-level paths in the reply.

"Top" means top of Finder's Date Added sort: `touch` does NOT change Date Added (2026-11-05 a re-shown map stayed buried and the
user said "its not there"); re-showing an old file needs a fresh copy under a new name (cp; check `mdls -name kMDItemDateAdded`).
