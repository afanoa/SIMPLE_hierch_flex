---
name: refine3d-shc-virgin-field-trap
description: "refine3D shc does the global greedy first pass only if the WHOLE ptcl3D field is virgin (every record e1=e2=e3=corr=0, NOT updatecnt); the 10180 327k \"realign\" (allptcl_realign, Sep 14) never searched globally and its poses are unconverged (median 112 deg off for the 46k good particles)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T18:48:15.660Z
---

**The trap.** In `simple_strategy3D_matcher.f90`, `has_been_searched = .not. is_virgin_field(oritype)`. `ori%has_been_searched` (simple_ori.f90:1039) is TRUE when e1, e2, e3 or corr is non-zero, NOT based on updatecnt. So a field is virgin only when every record, deselected rows included, has e1 = e2 = e3 = corr = 0. Zeroing updatecnt alone did nothing: the first retry ran 10% greedy and scanned 10% of the search space.
- So aligning new particles inside a project that already holds refined particles makes `shc` skip the greedy global pass for everyone. It starts a stochastic local search from whatever poses the new particles carry.
- `refine=greedy` is not a valid refine3D choice. The UI offers snhc, shc, neigh, shc_neigh, prob, prob_state, prob_neigh, pose_cont.

**CORRECTION (later on 2026-10-04):**
- `refine=greedy` IS honoured, although it is missing from the UI choices. The Sep 14 log shows mode greedy with 100% of the search space scanned; the "% GREEDY SEARCHES 0" counter only counts greedy steps inside shc.
- With `shc`, the workers see a non-virgin field even after e1/e2/e3/corr/updatecnt are zeroed (probably the master seeds orientations first), so shc runs 10% greedy and scans 10% of the search space. Use `refine=greedy` for a from-scratch global search.
- The Sep 14 failure was therefore NOT the search. Its reference was the abinitio recvol, which was "AUTO-RESCALED OLD-CONVENTION"; the reference is the suspect.
- Always include a blind control: re-align a sample of already-refined particles from zeroed poses and check they land within a few degrees.

**The consequence.** `10180/allptcl_realign` (refine3D maxits=1, Sep 14):
- 0% greedy searches, orientation overlap 0.025.
- For the 46,783 particles that built the reference, its poses differ from their abinitio/consensus poses by a median 112°. Only 2% are within 5°, and no single global rotation explains it.
- The "327k loses the conformational signal" observation and the "12 of 20 reproducible components" result were both on these poses.

**How to apply:**
- Before aligning a subset, zero e1/e2/e3/corr (binoris columns 8/9/10/3) and updatecnt (25) for EVERY record; state is column 23 and x/y are 27/30. Restore the kept particles' poses from an exported table afterwards.
- Verify iteration 1 reports ~100% GREEDY SEARCHES.
- Tools: `10180/autostates_c46k/emb327k/{patch_states.py,merge_poses.py,align281k.sbatch}`.

Related: [[spliceosome-realign-deselected-unlocks-basis]], [[spliceosome-region-flex-2026-10-04]].

**Blind control result (2026-10-04 15:23):**
- Run: refine3D refine=greedy maxits=1 volrec=no lp=8 trs=10 mskdiam=400 sigma_est=global, vol1 = the 5.91 Å 46k consensus, master build f9bbf72ce.
- Control: 2000 of the 46k refined particles, re-aligned blind from zeroed e1/e2/e3/corr.
- They land a median 87° from their consensus poses (projection direction 61°, in-plane 32°); only 9.6% are within 5°.
- New-pose score 0.371 vs consensus-pose score 0.372.
- Mirror and in-plane conventions don't explain it.
- So one-pass blind greedy does NOT recover poses for this dataset, the same failure as Sep 14. The merged project is renamed `p327k_merged_INVALID_control_failed.simple`.
- Open question for the user: how should the 281k be aligned instead (abinitio3D-style schedule, refine3D_auto from scratch, or iterative prob)?

**Control diagnostics (2026-10-04 evening, `emb327k/ctrl_diag_{A..G}`, 2000 refined particles, blind):** median error vs consensus poses.

| run | search | band | directions | other | median error |
|---|---|---|---|---|---|
| A | greedy | 8 Å | 2500 | sigma2 bootstrapped | 87° |
| B | greedy | 6 Å | 5000 | consensus sigma2 | 90° |
| C | greedy | 6 Å | 20000 | | 89° |
| D | greedy | 10 Å | 20000 | | 89° |
| F | greedy | 20 Å | 5000 | | 79° |
| G | prob | 8 Å | 5000 | | 81° (18% within 5°) |
| E | neigh (local), started FROM consensus poses | 6 Å | 5000 | | 0.8° (96% within 5°) |

- The wrong blind poses often score HIGHER than the true poses (F 0.384, G 0.376 vs 0.372).
- So refine3D's per-particle score prefers wrong orientations globally, even at 20 Å, while the true pose is a stable local peak.
- Suspects, not yet tested: CTF or stack-normalisation mismatch (would flatten the global landscape but leave local peaks), a defect in the blind-start path, or the reference/mask.
- Blind alignment of the 281k with refine3D is NOT possible until this is understood.
