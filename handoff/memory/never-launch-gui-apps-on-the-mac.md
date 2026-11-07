---
name: never-launch-gui-apps-on-the-mac
description: "never launch Chimera or other GUI apps on the user's Mac to render pictures; render headless (batch2 matplotlib) or give the user the command"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-03T21:33:32.234Z
---

On 2026-10-03 I tried to render region-coloured maps by launching Chimera.app with `open -a` and with `--script`. The user stopped me: "stop opening stuff".

**Why:** launching GUI apps takes over the user's desktop.

**How to apply:**
- To make pictures, render headless on batch2 (Python/matplotlib in an sbatch job) and send the PNG with SendUserFile.
- Or give the user a one-line ChimeraX command to run themselves.
- Never `open` an app, and never start a GUI binary locally.
