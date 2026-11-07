---
name: half-maps-always-in-halves-subfolder
description: "every map folder the user may open must hold FULL maps only: even/odd half maps go into a halves/ subfolder (tidy_halves.sh), automatically at the end of every job; user had to ask many times (2026-10-20)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-20T19:27:27.154Z
---

2026-10-20 the user asked repeatedly ("I can't open all these even odd volumes at the same time. I have asked you many times
to sort out the even and odd volumes") - consolidation outputs write recvol_stateNN{,_even,_odd}.mrc interleaved into
fam0/regionall/level1/maps/ (states_dir has full maps only, but the user opens maps/ too).
**Why:** they open whole folders in a viewer; interleaved halves triple the files and make the folder unusable.
**How to apply:** every job that writes state maps ends with routes/toy55/tidy_halves.sh <consolidation dir> (moves
*_even/_odd into maps/halves/; landscape_sample, merge_up, mu_family already do); readers that need halves must look in
maps/halves/ too (landscape_sample.py, compose_state_maps.py do). When pointing the user to maps, give the folder that holds
only full maps, and check it before saying so. Act on the run/folder they name ([[act-on-the-literal-ask]]).
