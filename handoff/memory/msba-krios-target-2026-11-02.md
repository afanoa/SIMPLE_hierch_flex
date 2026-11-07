---
name: msba-krios-target-2026-11-02
description: "user 2026-11-02: MsbA (membrane protein) is the target for the tree + pose-instability question; their own Krios MsbA data is being picked on batch2 (jobs msba_krios_refpick / simple_msba_krios), so no EMPIAR download; TMEM16A stays the published-data fallback"
metadata:
  type: project
---

User, 2026-11-02, while the mic_het fixture verdict was in: "msba is a membrane protein that will be a good target". The squeue of
that hour shows their own MsbA Krios processing on batch2 (msba_krios_refpick, 10 x simple_msba_krios workers on csbdevel), so the
real membrane data for the pose-instability question (mic_het said: no selective statistic, the guard's locality decides, fixposes
only as good as the input poses) will be the user's MsbA set, not an EMPIAR download.
**How to apply:** when MsbA particles exist, run the tree on them with the mic_het readouts (pose statistic per node, guard L/L_ref,
state maps) before any synthetic follow-up; keep [[tmem16a-empiar10124-target-2026-11-01]] as the published-data fallback.
Related: [[micelle-hetero-fixture-mic-het-2026-11-02]], [[tree-pose-artefact-test-2026-11-02]].
