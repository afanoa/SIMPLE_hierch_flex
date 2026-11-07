---
name: took-off-site-class-is-no-pole-2026-11-02
description: FIX (cb3196f0a on claude/tree-lever-floor): a site_rank endpoint whose EM class grew past half the node is no pole and claims no seeds; it had taken the real pole's seeds at the 8 % detlim root (717 vs 59,145) and backed the nulls' quality axis; gates identical to HEAD to the particle, detection limit unchanged (16-17)
metadata:
  type: project
---
2026-11-02, session tree-lever-floor (note hetero_tree_reworks_2026_11.md 9.1/9.3; batch2 math_reworks/lever, arms head / bc):
- Cause of the seed anomaly: seed_poles gives a particle to the pole with the highest posterior >= 0.5; a took-off class
  (share 0.83-0.97, posterior ~0.9 everywhere) wins against a real pole's 0.6 nearly everywhere, so the rare pole got 717 / 545
  seeds and the pole rule rejected both. The halves test's 18,600-seed pole on the same fits was an accident of averaging.
- Fix: merge_poles_subspace takes l_ok = class share <= 0.5; a rejected rest joins its best pole only while it is a minority
  (a majority rest stays the complement, the took-off pole's old role); the pole rule ignores a rest without endpoints.
- Gates: nulls 1 leaf; r131 both seeds identical trees to HEAD (corner kept); OLE same root/leaves (95 particles at r3); qtoy2 root partition moves <= 200 particles
  (sites 3-5 lose took-off-backed second poles); 8 % root pole {z1lo,z3hi} 1,863 seeds at purity ~0.75 but seed z 2.2-2.3 < 3
  -> still the fallback; limit unchanged. Rejected forms: the complement (flips the r131 root pole) and the climb (see
  [[site-climb-lever-no-go-2026-11-02]]).
**How to apply:** took-off endpoints are a diagnostic of a start that found nothing; never re-admit them as poles or complements.
Related: [[rare-state-detection-limit-2026-10-29]], [[math-reworks-verdicts-2026-11-02]].
