---
name: peer-message-to-idle-session-can-sit-unread
description: "TRAP: a SendMessage (uds socket) to an IDLE peer session can land in its transcript without starting a turn; check list_events after about 30-60 min of silence and resend with ccd send_message, which reports 'delivered; its turn has started'"
metadata:
  node_type: memory
  type: reference
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-24T11:21:54.188Z
---

2026-10-24 about 05:35 I sent the CAHRA session (local_6a17d89d, "Flexible molecules for tree flex model") a request for three
readouts via SendMessage to its uds socket. The tool said "in that session's inbox, not yet read". No delivery notice followed,
and two hours later its transcript (list_events) ended on my message as a [user] turn with NO reply: the session was idle and
never ran a turn on it. The same session had answered earlier messages on the same route within minutes.

**Fix:** mcp__ccd_session_mgmt__send_message with the local_... session id returned "delivered; its turn has started on it".

**How to apply:**
- When a peer owes me a result and stays silent past its expected time, read its transcript tail with list_events (limit about
  8) before waiting longer.
- If my message is the last turn with no reply, resend through ccd send_message (by session id) and log it.
- A peer missing from ListAgents may still be alive under its session id.
- Related: [[relay-intermediate-results-from-background-schedulers]].
