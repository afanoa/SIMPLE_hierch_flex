---
name: say-upfront-when-i-cannot-do-the-ask
description: "when the user asks for something I am not allowed to do (permanent deletion / emptying the Trash), say so plainly in the FIRST reply with their one-step way to do it; never quietly substitute a weaker action"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1a1d7175-9257-486d-a2bb-2f28c38b4e51
  modified: 2026-10-07T14:52:55.534Z
---

2026-10-07: the user asked twice to delete downloaded map folders on the Mac. I moved them to the Trash
(permanent deletion is outside what I may do) and mentioned it only at the end of a status reply; the
user had to ask again and then: "so are you unable to empty the trash. i have asked multiple times idiot".

**Why:** a quiet substitute reads as ignoring the request and costs the user a round trip.
**How to apply:** first line of the reply: "I can't permanently delete files or empty the Trash; I moved
X to the Trash, empty it in Finder (Finder > Empty Trash)". Also: this shell cannot read ~/.Trash
(macOS privacy), so I cannot even verify what is in it. For downloads the user may want gone later,
mention the size and location so they can remove them in one step.
Related: [[act-on-the-literal-ask]].
