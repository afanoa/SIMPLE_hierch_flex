---
name: tree-must-be-input-free-2026-10-28
description: "STANDING (user 2026-10-28): hetero_tree must be input-free - never seed a minority state from prior knowledge; references (closed rows, truth) only score. Scientifically positive, input-free changes: implement without asking"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T18:56:18.462Z
---

2026-10-28 ~15:10, after I called route A "the closed-state seeding run": "well are approach needs to be input free, i.e. we cant
seed a minority state if we dont know its going to be there. If these are scientifically positive stuff, feel free to implement."

**Why:** the tree has to find states on new data where nobody knows the answer; anything tuned to a known state is not a method.

**How to apply:**
- Every seed, split and merge must come from the data at that node. The 8,720 closed reference rows, qtoy2 truth and null labels
  are for SCORING only. Route A is input-free: `state_seeds` (simple_hetero_tree.f90) seeds from every discovered site pole that is
  not rank-correlated with the seeding site (< FGATE), has >= minpop pure members and stays < 1/2 node. Never call it "closed-state
  seeding" - say "extra seeds from discovered site poles".
- Do not adopt a change that only helps because it was pointed at a known state (e.g. early test 0 that seeded the closed pole by hand).
- Implement input-free changes that are positive on the gates without asking ([[science-before-speed-2026-10-28]]); bring mixed
  ones (a gain on one gate, a loss on another) to the user.
