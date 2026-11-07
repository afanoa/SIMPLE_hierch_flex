---
name: max-two-comment-lines-per-change
description: "user 2026-10-16: \"ridiculous. No more than two comment lines\" - a code change adds at most two comment lines in total; no routine headers, rationale paragraphs or trailing !< notes"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 950a20f7-af3a-4db4-93f0-bde563ac38a8
  modified: 2026-10-16T16:58:31.306Z
---

2026-10-16, flex PCG memory change (simple_flex_pca_pcg.f90, +513 lines): I had added 54 comment lines + 3 trailing
`!<` comments (same ~10 % density as the file). User: "how many comments did you add..." then "ridiculous. No more than
two comment lines."

**Why:** the user reads the diffs; explanations belong in my reply or the commit message, not the code.

**How to apply:** a whole change adds at most two comment lines. No per-routine `!>` headers, no rationale
blocks, no trailing `!<` on new declarations. Fix an existing comment the change makes wrong (counts toward the two),
restore existing comments I reworded, and move existing comments with their code unchanged. Count before reporting:
`git diff -U0 <base> -- src | grep '^+' | grep -E '^\+\s*!' | grep -vE '^\+\s*!\$'`.
Related: [[act-on-the-literal-ask]].
