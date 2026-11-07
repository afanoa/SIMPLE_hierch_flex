---
name: choose-routing-fixes-closed-dilution-2026-10-28
description: "state_split children are the cross-half EM's REGROUPING (coarse <=20 A), not the stage-B states -> closed state s02 delivered fused with s01; 'choose' routing (stage-B states when they win by > 3 SE on the held-out band) delivers s02 intact, nulls unchanged; ADOPTED default-on in recipe14"
metadata:
  type: project
---

2026-10-28 ~18:00, adaptive split/merge session, 95 replays / 28 core-h (beegfs hetero_tree_20261028/route9/, route9_on_v5.diff,
design doc section 10.15 in its worktree):
- Cause: the re-test's cross-half EM regroups on <= 20 A detail and moves 33 % of r111's particles off their stage-B state; existing
  route modes map states to the majority group, so both states landed in one group.
- choose = stage-B states delivered when they beat the regrouped groups on the held-out band by > 3 SE (won by 3-60x the bar on
  RNase P/OLE/qtoy2). s02 intact at route A r111 (25,528), base fixture r1 (15,139), strongest r1 (15,814). spost/band_direct fail.
- Guards: nulls stay 1 group (membership never changes K*); qtoy2 root same, c18 at r1211 recall/purity .97/.91 -> .90/.98; m33 closed
  pair kept as its two stage-B states. Old-reference purity up (.087 -> .110) but 73-76 % of reference rows kept vs 91-93 %.
- Open: OLE cross-seed reproducibility (why stage-B routing was dropped 2026-10-25) - same-rows r1 replay on/off (ole-repro).
- User 2026-10-28 rule: input-free + gate-positive -> implement; recipe14 = recipe13 (recipe12b + route A) + choose default on.
Related: [[split-routing-em-vs-stageB-2026-10-23]], [[site-pole-closed-state-replay-2026-10-27]], [[tree-must-be-input-free-2026-10-28]].
- OLE r1 same-rows check (recipe14, seeds 1-4, 91.5 core-h; beegfs choose_r1/): split ARI mean (min) em 0.87 (0.76) vs choose 0.82
  (0.66); stage-B map corr 0.98 (0.95) shared. Seeds 2-4 (XD K=6) agree 0.98 in both; only seed 1 (K=3 near-tie, strongest pole
  site 3) drops 0.76 -> 0.66 - EM regrouping was pulling the odd seed toward the others. Delivered maps reproduce 0.95-1.00.
  Verdict: kept (maps unchanged; loss confined to the XD-K near-tie seed).
- IgM (armsigm1 200k, m69M tree, 11 split nodes replayed; beegfs igm_r14_vs_r12b/replay/): choose delivered stage-B at 8/11 (3 under
  3 SE -> EM, identical). Of 7 motion nodes: dominant single-arm eta^2 sharper at 5 (+0.013..+0.17; r2 resolves arm 4 from a mixed
  arm2+4 EM split; r122 a5_th .293 -> .342), equal at 1 (r111), weaker at 1 (r121 a3_ph .211 -> .173); other-arm leakage down at 6/7.
  Junk split r21: +2,104 junk (+28 %) still pure. Extremes (10 % tails at purity >= .5) 0/20 at single splits in both arms.
