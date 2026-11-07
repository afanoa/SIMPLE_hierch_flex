---
name: site-shift-profile-verdict-2026-11-07
description: "VERDICT: a FLAT per-particle in-plane shift profile in the site classifier (xstageb) kills the IgM core-patch shift split (z 26 -> 1-5, R2 .32 -> 0) and keeps a1 DISCRETE (NMI .146 -> .129); a 1 px prior fails because the site's own shift information (~1 px^-2) is no better than the prior"
metadata:
  node_type: memory
  type: project
  originSessionId: 05f56e75-ddd5-4446-a0bd-dac2a2e73d7d
  modified: 2026-11-07T18:37:20.113Z
---

2026-11-07 (side session of local_e18b808a). Write-up: `igm_cryolo/tree/fixture_igm/xsb_shiftprof/VERDICT.md` on batch2; the classifier
is `xsb_shiftprof/scripts/xstageb_sp.py` (a copy of xstageb.py, env SPMODE proj | grad, SPREF group | mean, SPSIG 0 | px).

- **Flat profile works.** Two shift columns per particle are Schur-complemented out of every group score.
  - Core controls: z 26.7 / 25.8 -> 4.7 / 1.4 (swap / indep); site P/N 16-18 -> 1.2-1.7; label R2 on the shift error .32 -> .000.
  - a1: z 38.1 -> 33.8; NMI .146 -> .129; GAPw .10 -> .04, DISCRETE.
- **Tangent form matters.**
  - Use the masked in-plane GRADIENT of the groups' MEAN reference: T = sum_j R[j,u] P(M d_j mu). The mask stays put, as in
    M (R_k - mu).
  - The literal d/dx P(M R_k) per group also kills the core split, but costs a1 more (NMI .114). It has a mask-edge term (cos .86),
    and its per-group tangents differ by the displaced site itself.
- **A ~1 px prior fails** (core z 12, R2 .14-.17). The site's shift information G is ~.5-1 px^-2 at the core patch and .03 at an arm,
  so the site cannot measure a particle's shift better than the prior knows it.
- **Mechanism.** 93 % of the core split's group difference lies in the site's own shift span; 9 % of a1's does.

**Why:** settles the "profile the site's shift derivatives" lesson from [[pose-tangent-deflation-verdict-2026-11-07]].
**How to apply:**
- If the site path gets shift hygiene, use a flat profile of the mean-reference masked gradient. Never a prior at sub-pixel scale.
- A swap-design cross-half z of 2-5 on noise-level splits is a floor: swap couples the halves' noise, indep stays ~0.
- Untested: a2-a5, the null arms, a Fortran port.
Related: [[igm-fixture-xval-stageb-gapw-2026-11-06]], [[site-split-port-traps-2026-11-07]].
