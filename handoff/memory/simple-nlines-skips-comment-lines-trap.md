---
name: simple-nlines-skips-comment-lines-trap
description: "TRAP: SIMPLE's nlines() already skips '#'/'!'/';' lines, so nlines(f) - 1 for a header drops the last data row"
metadata:
  node_type: memory
  type: project
  originSessionId: d63181a2-43e0-44aa-a87c-6f7c55606053
  modified: 2026-11-08T07:14:57.861Z
---

`nlines` (`simple_fileio`) counts only lines that `str_is_comment` does not flag, and blank or `#` / `!` / `;` lines are not
counted. A table with a `# header` line therefore has `nlines = rows`, not `rows + 1`. Reading `nlines(f) - 1` rows after the
header silently drops the LAST particle.

It happened twice: `cheap_children` (2026-10-21, see [[state-tree-mergeup-port-2026-10-21]]), and the site chain's `chain_blocks`
coords read (v5-v8, 2026-11-08). The dropped particle got seed 0, i.e. no initial group in `cross_em`; the fix was
`site_chain_coords_read_fix.patch`.

**Why:** the bug is invisible: one row in 50-100k changes no score, and no read error fires, because the loop simply stops early.

**How to apply:** when reading a SIMPLE text table with a header, use `nlines(f)` as the data-row count, or read until EOF. Grep
new code for `nlines(` followed by `- 1`.
