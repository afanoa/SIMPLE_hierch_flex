---
name: batch2-ssh-banner-filter-trap
description: "TRAP: stripping the batch2 login banner with grep -v '^  ' also drops every output line that starts with spaces (uniq -c counts, indented scorer lines); echo a marker and sed -n '/MARKER/,$p' instead"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T15:43:07.840Z
---

2026-10-23: `ssh batch2 '...' | grep -v "^\*\|^  "` (meant to drop the HHS banner) silently removed `sort | uniq -c` counts
("  30000 0") and every indented line of the junk scorers, so a seeds file looked empty and per-type results looked absent.
**Why:** the banner lines start with "* " or two spaces; so does a lot of real output.
**How to apply:** start the remote command with `echo ===START===;` and pipe through `sed -n '/===START===/,$p'`; never
filter the banner by leading whitespace. Related: [[small-fast-cluster-steps]].
