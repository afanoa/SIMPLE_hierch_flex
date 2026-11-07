---
name: show-cavgs-means-copy-the-mrc
description: "show me the class averages" = scp the latest cavgs_iterNNN.mrc to the top of ~/Downloads immediately; no montages, comparisons or analysis first
metadata:
  type: feedback
---

When the user asks to see class averages, copy the newest `cavgs_iterNNN.mrc` stack straight to `~/Downloads/` (top level) in ONE command and say the path. No PNG montages, cropping, or reprojection comparisons first.

**Why:** 2026-11-03, IgM solve2D. I made montages and then started a class-vs-reprojection correlation job. The user said "just show me the classes man. I am better than your comparison" and "3 minutes and you cant even show me the mrc. that is ridiculous".

**How to apply:** `ls -t .../cavgs_iter*.mrc | head -1`, scp it to `~/Downloads/<dataset>_cavgs_iterNNN.mrc`, then reply with the path. Offer analysis only if asked. Related: [[maps-for-user-go-at-top-of-downloads]], [[do-the-requested-job-not-diagnostic-sprawl]].
