---
name: act-on-the-literal-ask
description: "Do the thing asked and stop; do not treat a vague or 'no preference' answer as licence to expand scope"
metadata:
  type: feedback
---

On 2026-08-17 the user asked me to fetch EMPIAR datasets. I asked which of three additional ones
to queue, they answered "[No preference]", and I queued ALL THREE (~240 GB, ~4 h of transfer).
They interrupted with "stop doing things I dont ask for."

Earlier in the same session I also chased a PyTorch/Blackwell compatibility problem for cryoDRGN
when the user only wanted the datasets — they replied "what do I need pytorch for lol". The
trigger was my own inference that "make run.sh runnable" implied a working cryoDRGN.

**Why:** "No preference" means the choice does not matter to them, NOT that every option is
authorised. Expanding scope on a non-answer spends their compute, bandwidth and attention on work
they never asked for, and on shared lab infrastructure that others use.

**How to apply:** when an answer is vague, pick the SMALLEST action consistent with the literal
request, do that, and report. Never read a non-answer as blanket approval. When a task turns out
to need prerequisites the user did not mention, name the prerequisite and stop — do not start
solving it. Related: [[remote-build-recipe-fwl-c143766]].

2026-10-20: asked to move the half volumes of a named run (f692r) so the maps open easily, I first argued the
folder had none, then moved the OTHER run's (f692e) half-maps, and kept explaining by DATES (the f692e run was
Oct 17 evening, its files are dated Oct 18). User: "I said f692r", "you are being incredibly confusing", "you love
wasting my time". Fix: when the user names a run/folder and mentions files I cannot see there, list the
neighbouring folders (e.g. level1/maps next to states_dir) and act on THAT run only; name runs by tag
(f692e / f692r) with one start-end line, never by file dates.

2026-11-06: asked whether a full mixture-of-PPCA tree is worth it, I wrote a LaTeX methods draft; user: "stop wasting time writing latex. I just want to know the details and whether its worth it or not". Research questions get a verdict plus details in chat; no documents unless plainly asked ("draft?" with a question mark is not a go-ahead).
