# Site layer -> tree states (session: site-layer assessment/integration, 2026-10-27)

User: skip the assessment, outline the integration. Outline: /Users/atifao/SIMPLE/doc/implementation_notes/site_layer_tree_integration_plan_2026_10_27.md
Nothing built, run or committed.

What touches other sessions (nothing applied yet):
- Route A changes advance() (the site step) -> seed_node in simple_hetero_tree.f90: all of a node's sites are ranked
  before seeding; a qualifying minority pole of a non-driving site becomes an extra seed state (K+m, not K x m).
- adaptive split/merge: the K* test should also judge an extra state in its own site mask; promoted / extra-seeded states
  should be fused at merge-up only on shown equivalence (or be exempt until that rule lands).
- port/perf: route A makes site ranks part of every tree, so the row-local site layer (your hand-over section 4b) is a
  prerequisite. The proposed base is d3513ff70, not the integration worktree's design layer.
- start-point: state_anneal would also apply to the extra state.
