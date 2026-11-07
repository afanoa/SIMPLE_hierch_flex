---
name: port-tree-site-layer-recovers-closed-2026-10-22
description: "factorial per-site layer on the port's full-depth simsites tree (fulltoy_cdisc, which misses the 8 % closed arm): ALL port candidates as sites recover closed partly (AUC 0.79 -> 0.86 after one truth-free re-mask; ceiling 0.95); bottleneck = discovery returns FRAGMENTS of the site mask; round-2 re-mask collapses the 2-class EM"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-23T09:06:02.619Z
---

Run: `port_tree_sites.sbatch` (tree_sites.py on a read-only view of a port hetero_tree, the port's own per-node candidates as the
discovery), `score_port_sites.py`, `cunion_em*.sbatch`, `site_mask_refine.sbatch` in worktree laughing-carson-4e8ca1, routes/toy55.
Outputs: rnasep_klem/factorial/port_sites_cdisc. Write-up: factorial_site_layer_results_2026_10_21.md, section "Site layer on the
port's full-depth toy trees".
- The port's discovery DOES produce closed-arm candidates; they rank below domain A by variance share and are fragments: 19-36 % of the
  true C each.
- As sites they pool to the root, independent of A. Best closed pole: precision 0.385 vs base 0.08, AUC 0.79. The tree's leaves sit
  at base.
- Mask tests (root site EM, 30 iterations, same seed):
  - true C mask: AUC 0.95 / top-2.5k 0.97;
  - union of the 4 fragments, chosen with truth: 0.91 / 0.91;
  - truth-free re-mask from the site EM's own two classes: round 1 0.86 / 0.81; round 2 (broader mask) the EM collapses into one class.
- Domain A as a site: precision 0.961 / AUC 0.92, better than the tree's leaf split. Junk sits in site poles at its base share; the tree
  removes it.

**How to apply:** keep ALL site candidates as per-site factors. The open problem is assembling a site's mask from its fragments, with a
guard against the 2-class EM collapsing.
Related: [[factorial-site-layer-test-2026-10-21]], [[per-node-junk-rule-prototype-2026-10-22]], [[pole-census-synthetic-validation-2026-10-22]].

**ADOPTED, LEAN SCOPE (user 2026-10-22 ~22:10).** The full scope (site_rank on all 188 full692 candidates, about 260 core-h at 32
threads, the root's 7 about half) drew "this sounds rather expensive...". The port builds the lean version into hetero_tree, on by
default:
- the tree is unchanged (site_01 drives the split);
- each distinct site runs once, at the highest node where it appears; Dice >= 0.5 with an ancestor's candidate = a copy, skipped;
- the extra site EMs fit on at most 100k of the node, then score everyone once;
- parent-level pooling EMs run only for non-copies;
- one re-mask round, only for pooled variables that pass the halves test, kept only while the pole class stays < 1/2 of the node;
- factorial output (label + posterior per site per particle).
Validate on fulltoy_cdisc, then run as a resume pass on full692. See [[price-new-stages-lead-with-lean-scope]].
- Evidence for the spec (factorial_site_layer_results_2026_10_21.md):
  - fit on a subsample, score all (`minor_em3s.py`): >= ~50k fitted gives Spearman >= 0.98 and top-N >= 0.90 vs the full fit -> the
    100k cap is safe;
  - the hold-out peak SE is 0.05-0.09 (`holdout_se.py`, NOT 0.01-0.02), so the refinement guard = no minority <-> majority flip plus
    no hold-out drop > 1 SE;
  - copy test = |Spearman| >= 0.5 maximised over the pole pairing.
- GATE VARIANTS (user: "test it overnight"; `score_variants.py`): V0 strict loses the 8 % closed site; V3 refine-first recovers it
  (EPM pole 0.40 / 0.31, AUC 0.79) with V0's false sites, so full692 runs V3 (job 54110499).
- Labels:
  - label = top-EPM (equivalent pure members), not top-N* (N* over-sizes diffuse poles to 101k);
  - EPM under-labels majority poles, so the majority = the complement of the minority WITHIN junk-screened rows (over all rows it
    falls to the base rate).
- OVERNIGHT RESULTS (2026-10-23):
  - full692 V3 finds the closed-vs-open site at r1 (r1_s5): top-2.5k maps read 0.84 / 0.52 (closed pole) and 0.49 / 0.99 (open pole).
    But the closed pole's labels are a swollen 148k complement (hard share 0.59), indistinguishable from A's genuine 70 % side, so use
    its posterior ranking.
  - IgM: all 5 arms tracked (|rho| 0.54-0.72); V3 = V0; 2 false sites.
