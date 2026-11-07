---
name: no-empiar-downloads-synthetic-igm-pdh
description: do not download the full IgM (EMPIAR-11077) or PDH (EMPIAR-10489) datasets; work on the synthetic arms/tethered stand-ins, the user will process the real ones themselves
metadata:
  type: feedback
---

User (2026-10-17): "dont download the full datasets. we can work on synthetic versions of those; I'll process them when I get a chance."
The challenge targets are IgM (five hinged Fab arms, EMPIAR-11077) and PDH (tethered E1/E3 at partial occupancy, EMPIAR-10489);
continuum sets like 10180/10028/TRPV1 are easy for flex_pca and only serve as controls.

**Why:** large downloads and processing are the user's call and time; the synthetic benchmarks
([[klem-synthetic-benchmarks-arms-tethered]]) carry truth files and are the development vehicle.
**How to apply:** develop and judge klem on `simarms`/`simtethered`; never start an EMPIAR fetch; when the user brings the
real IgM/PDH particles, run the chain on their project.
