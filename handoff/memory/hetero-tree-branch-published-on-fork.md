---
name: hetero-tree-branch-published-on-fork
description: "UPDATE 2026-11-07 (cleanup, [[fork-holds-two-branches-only]]): fork = hetero-tree e1d5fdd4c (gated, known to work) + hetero-tree-wip (ports + testing guide + handoff/); every other ref archived in the batch2 bundle ([[mac-retired-handoff-2026-11-07]]). Before: e1d5fdd4c (two docs-only notes commits on ee26d9bbf); before that ee26d9bbf (master merged into the tree, gated: [[tree-master-merge-gates-2026-11-06]]), fast-forward 72e5aaba0..ee26d9bbf (138 commits) at the user's word \"I need to push (to the fork today)\"; ports go on top. Earlier: UPDATE 2026-11-02: fork/hetero-tree = 72e5aaba0 (the combined code: step-1 cleanup, subspace site test, SVD discovery, draw-free split test, took-off pole, effect floor, numeric merge-up + its port, pose lines, notes and figures), fast-forward 9005a8c77..72e5aaba0 at the user's request; earlier: "2026-10-24: the hetero_tree branch tip d3513ff70 (worktree sharp-proskuriakova-b68280) is PUSHED to the user's fork (remote fork = afanoa/SIMPLE_hierch_flex, branch hetero-tree), so never rebase or amend at or below it; new work goes on top"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T23:07:59.748Z
---

2026-10-24, at the user's explicit request, the coordinator session pushed d3513ff70 to remote "fork" (git@github.com:afanoa/SIMPLE_hierch_flex.git), branch hetero-tree.
- d3513ff70 = the OLE code 16e52fbdd plus the population guard. The push was a fast-forward, 8421e1a31..d3513ff70. origin (hael/SIMPLE) was untouched.
- 2026-10-24 later, local only (not pushed), on top of it:
  - a027c6d03, lever 1: flex delivers its state maps distributed;
  - 8c1b66f31, lever (c): mergeup accumulate_jobs.
  Snapshot m8c62t = HEAD; m8c63t = HEAD + the default-off switches.
- Still uncommitted in the worktree (the user's pending decisions):
  - re-entry SIMPLE_TREE_REENTER (4b), now guarded once per lineage (a reentered.txt marker on the re-entered node);
  - the SPLIT_CONTRAST / ROUTE switches;
  - node_motion=yes (S_MOTION 15, motion_stats in simple_hetero_tree; the factor hand-back: a switch = two seeded poles A / B,
    a swing = A / middle / B, product cells folded under minpop), strongest_pole=yes and two_poles=yes;
  - the site_halves.txt -999.000 sentinel, and the site_pole_rule.txt post-rule state columns.
  Snapshot m8c69t (build 69, gate and 4 suites pass) holds everything; see [[ole-seed-split-choice-2026-10-25]]
  and [[package-tree-tests-2026-10-25]].
  - Later in the working tree (snapshots m8c71t to m8c75t, each gate and 4 suites pass):
    - the m8c71t seeding rule: a cell's core intersection when it holds minpop, else the whole cell;
    - the split-half gate (motion_gate, motion_floor, the matched-fake floor at |rho| >= 0.9 with permutation fallback);
    - motion_inherit, and motion_multi (no = one factor per node);
    - the coordinator's flex_pca component-sort patch (src/main/flex, applied with git apply from
      /Users/atifao/SIMPLE/.claude/patches/flex_sort_components_on_8c1b66f31.patch).
    - m8c77t: pole_rule fix (A) and the symmetry-mate fold. The Gamma patch is in the m8c78t snapshot only, not the worktree.
    - m8c79t (build 79, gate and 4 suites pass) = m8c77t + three changes in simple_hetero_tree.f90, all coordinator-approved:
      - the restart-path guard: seed_site(); S_RECB0 ranks site_<seed_site>, not site_01;
      - mirrored B: two_pole_seeds gains minpop. With no opposite end holding minpop, B = the |A| lowest kept-pole-score rows. The
        opposite search is skipped when the pole holds both ends of an axis. Hand-back lines log "B: <rule>";
      - empty fakes: no pole in either half = a null value of 0 (flag e). The null and floor logs count kept / dropped / empty / failed.
      The worktree and ~/SIMPLE_m8c hold exactly this. No tree runs on it yet; the m33 root replay is the test.
    - m8c80t = m8c79t + junk_screen forwarded to state_split and the merge-up (report-only under 'no'; both UIs register it).
      The m33 pair runs on it (full692_design80_s1, full692_base80_s1).
    - m8c81t = m8c80t + a developer split_ntests override on hetero_tree (replaces 2 N_root / minpop - 1). It is for the 4e attribution
      runs on cdisc (fulltoy_m81Mnt1_cd1 / m81Mnp0_cd1 / m81Mnt1_cd2).
    - m8c82t = m8c81t + SIMPLE_SPLIT_ROUTE=band (prototype env switch; frozen, not to be built into anything until the scorecard says routing matters).
    - m8c83t (build 83, gate and 4 suites pass) = m8c82t + speed only, results identical: (3) the site batch (S_SBATCH 18: under
      node_motion/strongest_pole every site of a node is ranked at once, sites/batch.txt); (9) under motion_multi=no the strongest
      candidate goes first, the rest run off the path after the hand-back (motion/offpath.txt, finished by off_work/S_OFFP; children
      wait for ancestors' off-path calls via anc_wait because inherit_cands reads ALL non-failing switch/swing of an ancestor).
      VERIFIED 2026-10-25 (cdisc M design maxdepth 2, fulltoy_m82Mdes_d2 vs m83Mdes_d2): identical sites / cands / seeds / split
      children / every candidate's call, verdict, bins and gate; only motion/offpath.txt extra; wall 3:16 -> 3:06 at nodepar 2.
    - m8c84t (build 84, gate and 4 suites pass) = m8c83t + SIMPLE_SPLIT_ROUTE=band_direct (prototype: when stage-B wins in the
      re-test band its states route whole, no full-band pass). Battery arm split_replay/m84_bdir (54279954). Frozen like 'band'.
    - m8c85t (build 85, gate and 4 suites pass) = m8c84t + SIMPLE_SPLIT_ROUTE=choose (prototype: one cross-half E-step on each RAW
      label set, per-particle re-test-band deltas; stage-B raw labels only if the margin > 3 SE (SE = sqrt(N) sd), else EM raw).
      cross_em gained niter/lpa/lpb, estep lpp. Battery arm split_replay/m85_choose (54283382).
    - m8c86t (build 86, gate and 4 suites pass) = m8c85t + motion_quality=yes (prototype key, default off; doc spec:
      laughing-carson-4e8ca1 doc/implementation_notes/tree_scorecard_2026_10_25.md "Self-calibrating cut"): node reference
      motion/qref (top vs bottom corr quintiles inside dfx deciles from recB0/b.simple, 2-state reconstruct3D, the bins' box),
      per candidate blur_db of the lower-mean-corr end bin vs the other from the *_unfil halves in the stable region; quality when
      dB >= cut/2 and >= 3 SE -> motion_call.txt 'quality (was ...)', neither chosen nor inherited. Replay script
      stree_port/stree_qual_replay.sbatch (copies a root WITH its bins). HELD 2026-10-25 15:00 until the factorial session confirms the
      pose/corr frame (stage-B poses gave dB_ref +177, n.simple poses -89 at the m33 root; the tree uses recB0 = stA0 frame).
      NOT BUILT further: the factorial validation failed (no rule separates quality from real candidates in node_motion's frame);
      the code was taken out of the worktree and kept as /Users/atifao/SIMPLE/.claude/patches/motion_quality_on_m8c85t.patch.
    - m8c87t (build 87, gate and 4 suites pass) = m8c85t + the PERMUTATION-NULL gate floor: null_work scores NULL_PERM=4 random-fifth
      splits per floor node (r, r1, r2) and nothing else; the fake-direction refits and their helpers (null_draw, write_fake_coords,
      real_posts, match_rho, no_pole, count_keep, append_file, start_names, tok_end, after_tok, NULL_DRAWS, GATE_MATCH) are gone.
      OLE floors under it: 0.025 / 0.020 (stree_port/floorrep), vs the fakes' 0.269 / 0.742. Coordinator-approved 2026-10-25.
    - m8c88t / m8c89t (builds 88-89, gate and 4 suites pass) = m8c87t + pooled_poles=yes (key, default off; site_rank's pooled_rule
      re-tests a merged group the member rule rejects on its members' AVERAGED half/quarter maps; merge_poles_halves gained lpool;
      the tree forwards the key to all three site_rank calls). m8c88t capped the margin at self when every partner was negative;
      m8c89t fixes it (self - best partner, as halves_accept_margin). Replays: stree_port/pooledrep (site_rank: stree_siterank_rep;
      nodes from recB0: stree_recb0rep2; one site's call: stree_s4call_rep). Coordinator rule: keep it as a candidate key only if
      the OLE r1/r2 cross-seed ARI (pooled vs control) holds or improves; never relax the share or seed-map tests to rescue r111.
      RESULT: holds. OLE r1/r2 ARI 0.714 / 0.735 for both control (m8c88t) and pooled (m8c89t), same children; flexsort m8c75t
      was 0.668 / 0.772. r111: both c2 ends pooled-accepted but the pole rule rejects (seed-map col 2.11, share 0.564). m33 root:
      site 4 becomes a swing candidate (margin 6.97, cos -0.279), site 6 still chosen. DECISION (coordinator): neutral in all
      three tests, so it stays in the code OFF by default and UNUSED (a change must improve the scorecard beyond noise; node-ARI
      build/seed noise is about +-0.05: the control alone moved 0.668/0.772 -> 0.714/0.735 from m8c75t to m8c88t). Revisit only if
      qtoy shows missed poles.
    - m8c90t (build 90, gate and 4 suites pass, 2026-10-25 18:56) = m8c89t with the flex sort patch REVERSED (git apply -R; src/main/flex
      equals HEAD 8c1b66f31 again: the original Gram / energy order; the sorted files are backed up in the scratchpad's
      pkg/flex_sorted_m8c89t/). THE BUILD FOR EVERYTHING from now on (coordinator): the sort lost the RNase P closed state at base80's root,
      so every RNase P result on m8c75t-89t is suspect ([[flex-component-sort-moves-site-seeding]]). The Gamma patch stays out.
- 4e, the Bonferroni re-test, is now a DEFAULT in the worktree (2026-10-24), uncommitted until the user's go:
  - test(ntests=): 50 shuffles, re-test threshold mean + z_crit * sd at alpha 1/21 / m, the board on the max of its first 20 shuffles;
  - the new state_split key split_ntests; the tree computes m = 2 N_root / minpop - 1 from node_r/rows.txt;
  - the merge-up calls test() without ntests, so it is unchanged;
  - the SIMPLE_SPLIT_NSHUF prototype switch is removed.
  Snapshot m8c64t = HEAD + the switches + 4e; gate and 4 suites pass.

**Why:** published history can't be rewritten without breaking the fork.
**How to apply:**
- Never rebase, amend or reset at or below d3513ff70 on this branch. New commits go on top.
- Pushing is still only on the user's explicit word ([[never-push-to-origin-without-permission]]).
- STANDING (user, 2026-10-24): "we push only to the fork". Pushes go only to the fork, never origin, after the batch2 gate.
  The coordinator does them; my commits stay local.
- 2026-10-29 CONSOLIDATED (user ask, local only, NOT pushed): local branch `hetero-tree` (worktree sharp-proskuriakova-b68280,
  upstream fork/hetero-tree, 14 ahead). On 8c1b66f31: fbfdb19bd build-90 working state (the old uncommitted worktree),
  6414aa450 recipe9 base, bbc2201a1 recipe12 levers, e3d02bca7 equiv merge-up default, 923b80a9b strongest_pole default,
  b2cc87b6e route A (13), 9aaf70855 choose (14), 393b451d9 locality guard (16), 50db57bdb kway_maxdepth (17), 0a22a8a5c
  row-local (17b), f439b4fe4 reuse/prune (17c), 5709bc083 STATES_DELIVERED (17d). Tip == beegfs recipe17d/source (rsync -c;
  only doc/code_overview + an .orig differ). Not built from git; the recipe17d build is the gated one.
- PUSHED 2026-10-29 (user: "push the latest commits to the fork"): fork/hetero-tree d3513ff70..5709bc083 (fast-forward, 14 commits,
  = recipe17d source; build recipe17d/build_r17d fast gate 15/15, 2026-10-28). Never rewrite at or below 5709bc083 now.
  Production run settings are env, not code: SIMPLE_SPLIT_CONTRAST=maps, SIMPLE_SPLIT_RHOMIN unset ([[floor-off-guard-verdict-2026-10-29]]).
- 2026-10-30 (user: commit the guard fix, build, fast gate, don't push): local hetero-tree = 21fac34db (locality guard children path,
  locality_guard_children_r17d.diff d6e45d9f) on 5709bc083, 1 ahead of fork. Built from a batch2 worktree checkout via SLURM job
  54722780 on norm (batch2 login node was down): hetero_tree_20261028/recipe17d_guard/build_guard2, fast gate 15/15. Not pushed.
  The aborted login-node build left a partial recipe17d_guard/build_guard (safe to clear).
- PUSHED 2026-10-30 (user ask): fork/hetero-tree 5709bc083..21fac34db (fast-forward).
- PUSHED 2026-11-01 (user ask, after build12 fast gate): fork/hetero-tree 21fac34db..9005a8c77 (tree distribution, site_em, junk cut; branch claude/tree-dist). Never rewrite at or below 9005a8c77.

UPDATE 2026-11-02 19:25 EST, user: "merge the lever-floor branch and the merge-up port into claude/tree-dist and push it to the fork".
claude/tree-dist was rebased (three note commits) onto the merge-up port head 55ee9c815 (= lever-floor c5176bddc + c569be829 + 55ee9c815),
source byte-identical to build_mu's (fast gate 15/15, seven gates, 692k merge-up re-run) and, minus the port, to build_lf's (the 692k run
f692_lf_s1). Pushed fork hetero-tree 9005a8c77..72e5aaba0 (27 commits); origin untouched. NOT on the fork: the core budget
(claude/objective-meitner, dc72c27fc, default off), the prototype notes (posecond, ppinfo), the mic_het fixture branch.
