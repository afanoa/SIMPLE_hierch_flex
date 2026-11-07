---
name: at-node-switch-swing-call-2026-10-24
description: "At-node switch/swing call calibrated (phase E, 52 series, 25.5 core-h): nb_SF >= 0.030 + residual z >= 5 = swing; switches 8/8, swings 22/23, but the tree's false splits are NOT map-null (composition gradients, small nodes ride on views/defocus); RNase P m33 all switch; spec in doc/implementation_notes/at_node_switch_swing_rule_2026_10_24.md"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-25T04:22:42.036Z
---

2026-10-24, for the coordinator. The user's direction: keep the tree, but at each split call switch (split it) or swing (pole states plus
an annotation). Tools are in routes/toy55: phaseE_maptests.py, phaseF_rel.py, view_check.py, calE_final.py, and the jobs calE_*.sbatch,
rnasep_m33.sbatch. Outputs are in rnasep_klem/factorial/phaseE/; results in factorial_site_layer_results_2026_10_21.md, phase E.
- Bins are cut at the node: 5 quantiles of the candidate's score over the node's rows, tree junk out, no truth. Scores:
  - kept pole: the mean of its endpoint posteriors;
  - fflex fallback: the Fisher discriminant of the node fflex z between the seed groups.
- The call is nb_SF, the reproducible non-blend residual over the reproducible end power in the focus mask.
  - The residual/noise R mixes size with particle count: R - 1 ~ nb x ends SNR. RNase P r1 reads R 18-22 at nb 0.011.
  - The cut is 0.030: switches 0.004-0.025, swings 0.035-0.52.
  - The one miss is m8c17 a5_th at r111 (4.6k / bin): a small-amplitude arm that reads as a blend.
- The false splits are reproducible local map differences. z_ends, local change size, cosine to the mean map and FSC all overlap the
  real changes.
  - cdisc's null-type nodes are junk / composition gradients that blend perfectly, so they call switch.
  - The small IgM nodes follow views or defocus (end-bin view TV excess 0.44-0.64; r111212 is a defocus split).
  - All-junk r2122 (811 / bin) is the one false swing.
- CAHRA is NOT brightness. Contrast vs open is -0.07, and P5 - P1 has cosine -0.14 to the mean map. The focus fraction rises with SNR,
  which corrects [[latent-graph-cannot-tell-discrete-from-continuous-2026-10-24]].
- RNase P m33: every motion at root, r1 and r2 calls switch (nb_SF <= 0.013), most of them graded. r2 has two separate motions.
- OLE reference for the port's node_motion (2026-10-25, ole_12707/motion_calls/, box_crop 200): all 9 r1 / r2 candidates over both
  seeds call switch (nb_SF <= 0.004). The shared axes agree across seeds. The r2 axis, cut at opposite ends by the two seeds, is an
  occupancy switch, not a swing. The graded flag flips across seeds near its 0.40 cut.
- TRAP: m8c17's kept sites are scored on their own levels, with no population-wide column (the root level covers 81 % of rows). Use
  pick_igm2.py with node-level projects.
- TRAP: when the last-ranked site gives no split, the seed dir can lack site_pole_rule.txt (the fflex fallback). null_score2.py handles
  it.
- TRAP (2026-10-25): site_pole_rule.txt numbers its rows in its own seed order, which can put the rest first (IgM M fsite_r1_s3 keeps row
  2 = pole line 1). Take the pole from the seeds the tree split on (site_seeds.txt = node seeds) and its site_poles.txt line.
  null_score3.py does this. All 52 series were checked and none is affected.
**Why:** the tree needs a per-candidate switch/swing decision. Calibrated this way, it separates motion types but not real from
spurious splits.
**How to apply:**
- Port the spec as written.
- Don't expect the map call to remove null-type splits; that needs the junk / composition model or a recurrence test.
- Report the acquisition flag (view excess >= 0.40, defocus excess >= 0.20) rather than gating on it; its margins are thin.
- Related: [[discrete-splits-capture-the-span-states-that-matter]], [[view-imbalance-inflates-between-state-metric]].
