---
name: contrast-switch-kills-real-rnasep-splits-2026-10-27
description: SIMPLE_SPLIT_CONTRAST=yes removes 3 real unique m33 RNase P splits (not build/chance); the null qtoy splits are tiny real differences at huge N (no leak); relative-effect rho separates nulls (<=0.023) from real splits (>=0.039) with contrast on
metadata:
  node_type: memory
  type: project
  originSessionId: 01619764-7985-49aa-b2e9-539d9cd5a16f
  modified: 2026-10-27T21:13:21.546Z
---

2026-10-27, adaptive split/merge session (test builds ~/SIMPLE_smr4 -> ~/SIMPLE_smr5o, login-node; outputs
ole_12707/split_merge_replay/{emnull,hpmask,ortho,rho,arms}/):
- m33 r1211112 / r2111112 / r2121212: K* 2 on m8c36 = d3513ff70 = test build with contrast unset (z 3.3-3.6, rho 0.11-0.25);
  K* 1 only with SIMPLE_SPLIT_CONTRAST=yes (z -0.45/1.24/1.20). Deterministic. The parent showed the 3 differences are real,
  reproducible (eo 0.65-0.70) and carried by no other state pair. So contrast=yes costs real RNase P states.
- qtoy null featureless root split (z 5.2 with contrast=yes): no leak (cross_em is truly cross-half; high-pass before mask and
  per-shell orthogonal difference leave T unchanged); EM-inclusive null is degenerate (random starts drain). It is a tiny real
  difference at ~30k/group (misaligned-particle share 3.7 vs 6.3 %).
- rho = difference SNR / consensus SNR per particle in the re-test band: with contrast=yes nulls 0.012-0.023, real/positive
  splits 0.039-1.47 (m33 closed pair 0.046). Without contrast the qtoy null quality root is ~0.1 (estimate) -> floor alone
  may not suffice; contrast-unset rho battery was running at session end.
Related: [[split-test-accepts-quality-split-on-null-2026-10-27]], [[rnasep-closed-state-must-be-retained]],
[[equivalence-mergeup-verdict-2026-10-27]].
