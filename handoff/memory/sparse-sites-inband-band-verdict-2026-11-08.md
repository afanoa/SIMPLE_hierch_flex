---
name: sparse-sites-inband-band-verdict-2026-11-08
description: "VERDICT: no in-model sparse-sites band (per-site ARD prior, per-block observation band by held-out gain, node covariance-FSC weights) recovers IgM 10 % fix at the full band; seeds need the band from their FIRST iteration (ladder caps do this); ladder stays"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-08T18:20:30.316Z
---

sparse-ppca side session, 2026-11-08. The user asked for the band inside the model instead of the ladder
([[sparse-sites-band-ladder-verdict-2026-11-08]], [[band-must-be-data-derived]]).
- Note: `handoff/sparse-ppca/inband_band.md`.
- Patches and builds: `math_reworks/sparse_port/` (b8d, b9d, inb3, inb4).
- Runs: `native/ib{1,2,3,4}_*`.

**What failed, and why:**
- **v1, per-site per-shell ARD prior (MacKay ML-II) plus a prior-whitened group norm.** IgM: no site.
  - Every seed's learned edge stayed at the full band while its energy fell 14x, as without ARD.
  - The prior reads the mode's own per-shell energy. A weak seed's high shells shrink with its low shells, and the ball's hard edge keeps high-shell content, so the profile never cuts. The group penalty is not involved: phase 1 has none.
- **v2 / v3, per-block observation band.**
  - Each block's E-step, evidence and lasso use only shells up to its edge. The modes are fitted at every shell, and the edge is the argmax of the held-out gain per edge.
  - Box 36 is fine: v3 qtoyflex within .01, and r131 reads higher.
  - IgM seeds collapse faster than at the full band (32.7 to 0.003). Narrowing the E-step after one iteration starves the still-random seed's latent.
  - Traps:
    - per-edge LAPACK calls in a hot loop serialise in OpenBLAS (5-7x cost); batch them;
    - an all-negative gain curve sends the argmax to ring 2 forever; move the edge only on a positive gain.
- **v4, one node weight** (coordinator's design): the even/odd FSC of the node covariance along one power-iterated probe, turned into Wiener weights that scale the data.
  - IgM weights fall below .5 beyond ring 2-6. The single probe finds a global low-frequency direction, not the localized arms, and the 10 % fix finds no site.
  - Box 36: the down-weighting costs qtoyflex's best blocks .02-.05.

**How to apply:**
- Any in-model band must be fixed BEFORE phase 1, by a statistic that sees localized variability. One untested candidate: covariance even/odd agreement inside the seed balls, pooled.
- The dropped-round drift persists under every version.
- The ladder (band3) remains the working rule. It can be made cheaper without changing a pick: caps 1 and 2 as concurrent jobs, and a round-1 screen for box-36 nulls.
- Coordinator 2026-11-08: briefed the seed-ball even/odd band as the last in-model try. Whichever wins (it, or the ladder plus
  concurrent caps 1-2 and a round-1 null screen) is built, fast-gated and committed to wip.
