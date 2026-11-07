---
name: klem-oracle-labels-vs-board-numbering-trap
description: "TRAP: oracle_arm_init.py init_labels number cells g1 + n1*g2 + 1, the klem board's argmax uses m = i1*n2 + i2 + 1; comparing them directly gives ~chance even for a perfect board (the '5.7 % stay' oracle claim was this bug)"
metadata:
  type: feedback
---

The board's lattice index (klem_lattice_kernel): i1 = (m-1)/n2, i2 = mod(m-1, n2); cell files cell_<i1>_<i2>.mrc with
0-based i1/i2 = the "g1 g2" columns of init_positions.txt. oracle_arm_init.py's init_labels.txt use g1 + n1*g2 + 1.
**Why:** oracle_region.sbatch's "argmax cell == oracle start cell" compared the two and reported 4-11 % for boards
that actually kept 43 % exact / 85 % within one step. **How to apply:** compare in (g1, g2) coordinates
(oracle_agree.py), or report corr(z1, g1) / corr(z2, g2); never compare label integers across these files.
Related: [[igm-arms-ceiling-pose-limited-2026-10-19]].
