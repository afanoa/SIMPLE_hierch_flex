# hetero_tree paper draft + integrated-node prototype (2026-11-06, session magical-babbage)

User asked to start a paper (3DVA / RECOVAR positioning, what each molecule adds, new molecules) and then: "Put the integrated node
model in this paper: prototype it, then run the current tree as the baseline it is compared against. Start drafting the sections whose
results already stand."

## Draft
- LaTeX: worktree magical-babbage-701c7b, doc/manuscripts/hetero_tree_2026_11/main.tex (+ refs.bib, figures/, main.pdf; compile with
  ~/.local/share/micromamba/envs/tex/bin/tectonic -X compile main.tex). Same layout as the flex2D draft. Orange [..] = pending.
- Drafted: abstract, introduction, results 2.1 (overview), 2.2 (calibration), 2.3 (RNase P), 2.4 (OLE), 2.7 (limits), methods
  (node model, sites/poles, refinement, draw-free split test, guard, effect floor, integrated node maths, hinge fixtures, data sets).
  Pending: 2.5 integrated results (interim table in place), 2.6 3DVA/RECOVAR head-to-head, discussion prose, bib details marked [check].

## Prototype (batch2 /mnt/beegfs/elmlund/afan/paper_inode)
- scripts/inode_em*.py: mixture of PPCA per node, z integrated out (Woodbury), local pose search, rigid tangents projected off U_k,
  half-set designs swap / indep / merged (stage-B analog). v2 = --retest + nested held-out rival; v3 = --seeds_file, per-half refits in
  the retest; v4 = --contrast (per-particle contrast profiled on a grid). export_node.py turns a tree node into a toy-style folder
  (handles the RNase P stkind/indstk lookup and ctf=flip via --ctfabs 1).
- Toy hinge fixtures (toy3 39,990 ptcls; f50, f15 = 15 % minority, cont, null) + production-tree baselines (build_mu, root only) in
  trees/<fix>. Scores: logs/inode_score.54993121.out.

## First verdicts (seed 1, interim)
- Production tree root: f15 NMI 0.967 (2 leaves), f50 0.978, cont split at the median, null 1 leaf; root stage-B pose scatter 1.85 deg.
- Integrated node: f15 NMI 0.91-0.93, f50 0.95 (= stage-B analog, > fixed poses 0.86; true poses 1.0); pose scatter 1.48-1.68 (input 1.52,
  greedy stage B 1.77); within-state pose bias unchanged by every arm (~1.4 deg median).
- Held-out two-vs-one likelihood: prefers two states on every fixture incl. the continuum (1.2-2.0) and the null (0.17-0.70 nats/ptcl
  depending on the training design): not a calibrated split test yet.
- Next (submitted by a waiter when the export lands): RNase P node r1 of f692_lf_s1 with the tree's own seeds, integrated vs stage-B
  analog vs fixed poses, scored by the 8,720 closed reference rows (tree: 89 % of r1's closed rows into r13 at purity 0.063).
