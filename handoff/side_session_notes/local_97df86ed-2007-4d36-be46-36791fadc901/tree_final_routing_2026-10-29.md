# Final routing pass for hetero_tree (tree-final-routing session, 2026-10-29)

## What it is
- Code: worktree tree-final-routing-f34aad, src/main/strategies/search/simple_hetero_tree.f90 only (uncommitted).
  Diffs on the cluster: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/route/route_v1.diff (md5 aeebfcc9), route_v2.diff (ab92cc0a).
  Builds: route/build_v1 (gate 15/15), route/build_v2 (v1 + node scoping, pruning, flat iteration count).
- Env, default off: SIMPLE_TREE_ROUTE=only (offline on a finished tree: run() does only the route step) or =yes (end-of-tree step
  after STATES_DELIVERED). v2 adds SIMPLE_TREE_ROUTE_ROOT=<node id> (one subtree), SIMPLE_TREE_ROUTE_NODES=<id,id,..> (score only
  those routing nodes; flat skipped) or =flat (flat only), SIMPLE_TREE_ROUTE_FLATITS=n (default 2).
- Method: route/route.simple = the root project with every row's leaf-node stage-B pose, state = leaf (junk leaves included).
  Each problem (a routing node = node with >= 2 children, or flat<k>) reconstructs one map per child (tree membership; flat: per
  leaf / per flat<k-1> label) with reconstruct3D (box 88, halves), then scores each half's rows against the OTHER half's maps with
  stage B's refine3D (prob_neigh geom, lp 10, box 88, nspace 5000, 1 iteration, local poses; eo alternated inside the job so the
  end-of-iteration reconstruction has two halves). Routed child = argmax of prob_state_evidence.txt's per-state log-evidence;
  margin = best minus second (nats). Route-down passes rows to the child's routing node (rows.txt) or lands them (landed.txt).
- Outputs in route/: route_summary.txt, route_rows.txt (row tree_leaf down_leaf first_disagreement_depth min_margin, flat leaf and
  margin per iteration), down_labels.txt / flat<k>_labels.txt (same format as hetero_tree_labels.txt), node_<id>/scores.txt.
- Launcher: route/route_run.sbatch (mirror run dir with copied tables, nothing written into the tree). Scorers: route/route_qtoy.py
  (truth), route/route_closed.py (closed reference rows).

## qtoy2 (recipe17b gates_env/qtoy2_s1, floor on; route/q17b_v1): 10.1 core-h, 7 min, 25 GB (v1 does not prune; v2 does)
| labels | own leaf | NMI combo | leaf-majority right |
|---|---|---|---|
| tree | 1 | 0.706 | 0.835 |
| route-down | 0.829 | 0.717 | 0.841 |
| flat it1 | 0.796 | 0.770 | 0.878 |
| flat it2 | 0.781 | 0.789 | 0.893 |
- Best-leaf F1 per combo, tree -> flat2: c1 .54->.69, c2 .65->.75, c18 .87->.95, c22 .78->.91, c28 .94->.98, c43 .61->.83,
  c49 .81->.92, c101 .94->.98. Genuine rows flat2 moves (8,956): right .43 -> .77. Route-down moves (7,639): .53 -> .57.
- First disagreement depth: 0 6.1 %, 1 6.9 %, 2 3.7 %, 3-4 < 0.5 %. Median min-margin 5.7 nats agree vs 1.0 disagree; 16.6 % < 1 nat.
- Cost per step: flat iteration 2.5 core-h (60k x 13 maps); root 1.0 (60k x 4); ~0.3 core-h fixed per node.

## v2-v5 and the iteration curve
- v3 = v2 + SIMPLE_TREE_ROUTE_LEAVES=<leaf,..> (flat over a leaf pool). v4 fixes the node list (v3 scored every node when LEAVES was set).
  v5 = v4 + volrec=no on scoring jobs (refine3D's assignment-only switch: a state emptied in one half no longer kills the job with
  "carried-forward state volume is missing"). Gate 15/15 for v1-v5 (v3's gate rerun after a node failure). Current: route_v5.diff
  md5 55a61e66, route/build_v5.
- qtoy2 flat iterations (route/q17b_v5_flat4, flat only): NMI .770 / .789 / .809 / .827, majority .878 / .893 / .901 / .909, best-leaf
  F1 c1 .62/.69/.78/.82, c43 .73/.83/.89/.92, c2 .71/.75/.78/.80; own leaf .796/.781/.728/.716. 10.7 core-h, 13 min, 489 MB.
- Reproducibility (q17b_v1 vs q17b_v2_flat4): routed child identical >= 99.95 % at every node; flat1 99.97 %, flat2 99.76 %.
- Launcher: RESUME=1 reruns the controller in an existing mirror (finished steps are reused). The root project is taken from the node
  table (the integration tree has rnasep.simple, no input.simple).

## RNase P closed lineage (692k integration tree, route/f692_lineage_v3, v5): 34.7 core-h, 34 min, 7.1 GB
- Env: NODES=r,r1,r13,r131,r1311,r2,r23,r232; LEAVES = r13 + r23 subtree leaves (9, 114k rows); FLATITS=3.
- Own child: r .836 (refs .861), r1 .865, r13 .891, r131 .851, r1311 .788, r2 .752, r23 .846, r232 .814. Median margin ~2 nats;
  45.5 % of rows min margin < 1 nat. First disagreement: depth 0 16.4 %, 1 10.0 %, 2 1.4 %, 3 0.9 %, 4 0.4 %.
- Root, ref rows: r1 4447 stay / 307 -> r2; r2 2175 stay / 282 -> r1 (11 %); r1 keeps 4,987. r2-leaf refs in lineage leaves: 187 (down).
- r2322 refs (743): down keeps 434; -> r1222 155, r2321 70, r1312 39, r1311_s2 14, r1311_s1 8.
- Refs tree/down/flat1/flat2/flat3: r1311_s2 816/778/767/813/799; r1312 1824/1815/1664/1610/1551; r1311_s1 435/485/481/476/498;
  r2322 743/636/736/779/774. Lineage union recall .353/.353/.334/.332/.327, purity .069/.072/.066/.066/.065. The corner swaps ~40 %
  of its rows per flat iteration (about 5k in / 4.6k out) with ref counts unchanged.
- Arm readout (closed/open, scorecard recipe; rec1 = tree members in the routing recipe): r1311_s2 tree 1.00/0.58, rec1 1.02/0.60,
  after it1 0.97/0.59, after it2 0.97/0.59 (corner kept, 0.05 lower). r1312 .84/.56, .86/.58, .85/.60, .86/.60; r1311_s1 .77/.57,
  .78/.60, .76/.60, .76/.61; r2322 .78/.69, .80/.72, .74/.69, .75/.68. Distinct throughout.
- Margin gate (route_gate.py, flat1 scores): qtoy2 NMI .770/.763/.753/.741/.724/.709 at t = 0/1/2/3/5/10 nats (the gain needs the
  low-margin moves); RNase P lineage refs 2912 -> 3075 (tree) as t rises: the gate only removes the churn.
- Verdict: no default change; route-down = diagnostic; flat = optional for discrete-state data; full 692k flat (~300 core-h) not
  recommended on this evidence. Pending: floor-free qtoy2 KEEPMAPS tree with ROUTE=yes (route/qtoy2_keep_s1, job 54665026).
