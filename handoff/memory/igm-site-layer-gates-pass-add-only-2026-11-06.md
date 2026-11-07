---
name: igm-site-layer-gates-pass-add-only-2026-11-06
description: "per-site layer (islands + xval site stage B + split only on GAPw DISCRETE) replayed on saved nodes: 89 runs, 0 DISCRETE -> RNase P closed-state + null gates PASS but inert on RNase P; must stay ADD-ONLY (r21 / r211 per-site factors = the tree's own closed-lineage splits, GAPw calls them not discrete); 692k root and r5 untested"
metadata:
  node_type: memory
  type: project
  originSessionId: c79bc433-5d25-4389-9589-67d90cec0604
  modified: 2026-11-06T09:50:34.682Z
---

2026-11-06, IgM fixture session.
- Report: side-session-notes/local_5bd42cba.../tree_igm_fixture_report.md, section GATES.
- Work dir: batch fixture_igm/gates_x/. Scripts: gate_prep.sbatch, gate_xsb.sbatch, gate_submit{,2,3}.sh, gate_table.sh.
- GAPw calls: discrete_vs_continuous.md sections 7 / 7b.

Results:
- No gate run read DISCRETE, and no closed leaf was split into substates:
  - null fixtures: 24 runs;
  - RNase P r131 corner: 20 runs;
  - 692k closed lineage: 45 runs, over the leaves r1311 / r1312 / r1321 / r2112 and the ancestors r13 / r131 / r132 / r21 / r211.
- The strongest real-data per-site factors are the tree's OWN splits seen through one site:
  - r21 sites 1-3 match r211 | r212 (82-86 % agreement). GAPw: no call (G2 / G3).
  - r211 site 1 (r1 .97) matches r2111 | r2112 (83 %), the split that makes the closed leaf r2112. GAPw: not discrete (.45).
- The leaves pass partly for lack of power: GAPw needs >= ~50k particles even for the fixture's arm.

**Why:** the user asked for these gates before any port, and the RNase P closed state is a hard gate ([[rnasep-closed-state-must-be-retained]]).
**How to apply:**
- Keep the layer add-only. As the sole split rule, GAPw would stop at r21 and lose the closed leaf r2112.
- Untested where GAPw has power:
  - the 692k root (4 sites);
  - r5 (120k).
- No saved site masks (production prunes them) for r1322, r5111, r5121 or the r5 lineage.
- User 2026-11-06 skipped the Python ctun / ctpk / mm3q check and asked for the port proposal (notes folder: port_proposal_site_split.md).
  - Design: the layer acts only at would-be leaves above the depth cap; islands are for the layer only; the direct GAPw readout; site_split=no default.
  - Gates 0-1 come after approval; ctun / ctpk / mm3q move to Fortran gate 5. The ctpk and mm3q reruns have no stageB/b.simple.
  - f692 facts behind it: maxdepth 4; all 7 closed leaves sit at depth 4; 12 of the 66 leaves sit at depth 3, the only ones the layer can reach.
- Replay traps:
  - The island-gates reruns of r13 / r131 / r21 stopped before stage B. Use the production node's stage-B project with the rerun's sites.
  - The prep's node_rows.txt and the run's labels_final.txt are both in node order. Cross them with a child's rows.txt by a set lookup.
- Real IgM (683k consensus, ~1/3 of particles in the wrong 60-degree register): a DISCRETE call at the J slot or an arm site
  may be that register nuisance. Check it against register first ([[igm-11077-consensus-poses-tree-ready-2026-11-04]]).
- Related: [[igm-fixture-xval-stageb-gapw-2026-11-06]], [[discrete-vs-continuous-site-gap-2026-11-05]].
