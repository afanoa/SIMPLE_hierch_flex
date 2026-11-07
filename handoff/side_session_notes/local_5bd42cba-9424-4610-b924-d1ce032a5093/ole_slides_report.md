# OLE slides (2026-11-03)

Deck: https://claude.ai/artifact/TtsUmhA9q46jEiKfjvr7ub (10 slides, private until shared)

- Tree: 6_hetero_tree (seed 1, C2), 31 leaves (not 32). Leaf maps = final_maps/emp2 (~5 A).
- Internal nodes reconstructed fresh from consensus poses, C2, box_crop 256 (1.341 A/px):
  batch2 /mnt/beegfs/elmlund/afan/ole_12707/slides_20261102/recs/run01..11 (root 3.4 A, depth-1 3.6-3.9 A, depth-2 3.8-4.9 A).
  Display copies (low-passed at own FSC, solvent envelope): slides_20261102/maps/nodes, maps/leaves.
- Findings: cap of the handle 0.9-1.45 (contour units) in all r1 leaves, 0.1-0.3 in all r2 leaves, consensus 1.0.
  Lower C2 bridge helix absent in r11212, r222; weak in r1122, r1212, r221; consensus 1.2.
- cryoDRGN (cryodrgn_z8, existing run): clusters 21-68 % r2 particles, NMI to leaves 0.03; volumes carry the cap even in the r2-rich cluster.
  CAUTION: cryoDRGN volumes correlate only 0.65 with the consensus at 10 A in any axis frame (tree nodes 0.998), so no voxel-region numbers on them.
- Scripts: slides_20261102/{analyze2.py, render2.py, make_figs.py, nodes_prep.py, nodes_rec.sbatch}.

## Rework 2026-11-03 (user: plain grey, show motion, audience view)
- Deck rebuilt to 13 slides: molecule parts -> tree (split maps) -> cap -> bridge -> combination grid -> 4 branch slides (blink GIFs + split cards) -> every split real -> cryoDRGN -> appendix gallery.
- Split difference maps: slides_20261102/maps/splits (analyze3.py): child1-child2 from even/odd halves, 12 A, blobs >5 SD >=150 vox; all 30 splits pass, noise control 0 voxels everywhere.
- Handle rotation about the 2-fold (angles.py/angles.json): r1 +1..+3 deg, r2 -1..-4 deg vs consensus (~7 deg span); mostly a path/shape change, not a rigid spin.
- Overlays: png6 (shared grey, A-only blue, B-only orange); GIFs: gifs/.
