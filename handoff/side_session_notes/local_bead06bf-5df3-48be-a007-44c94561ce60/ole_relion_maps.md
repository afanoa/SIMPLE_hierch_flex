# OLE RELION maps for the deck (2026-11-05)

batch2 /mnt/beegfs/elmlund/afan/ole_12707/relion (README.md, states_table.txt, occ_table.txt).
All maps: box 360, 0.954 A/px, C2, FSC-weighted, B = 0 (unsharpened), unmasked, same frame as the SIMPLE tree maps (no mirror).

- 31 states (A = SIMPLE consensus poses + CTF-refined star, no local refine): PostProcess/cpose_w_<leaf>/postprocess.mrc
  = ~/Downloads/ole_relion_<leaf>_<res>A.mrc (halves Cpose/<leaf>/h1.mrc, h2.mrc = ~/Downloads/halves/ole_relion_<leaf>_half1|2.mrc)
- consensus (80,598, 3.15 A): PostProcess/cons_w/postprocess.mrc = ~/Downloads/ole_relion_consensus_3.15A.mrc
- r1 (r1* leaves, 50,606, 3.18 A): PostProcess/branch_w_r1/postprocess.mrc = ~/Downloads/ole_relion_r1_3.18A.mrc
- r2 (r2* leaves, 29,992, 3.40 A): PostProcess/branch_w_r2/postprocess.mrc = ~/Downloads/ole_relion_r2_3.40A.mrc
- B (local refine, worse in every state): PostProcess/leafd_w_<leaf>/postprocess.mrc = ~/Downloads/ole_relion_local_refine/
- Not for display: Cpose/<leaf>/full.mrc and ~/Downloads/ole_relion_unweighted_raw/ (unweighted to Nyquist, look like noise).

Slide regions (occ_table.txt, 8 A / equal-volume contour): cap and bridge track the tree maps across states (r = 0.99), at about half
the tree maps' spread. Cap in A: r1 states 2.1-2.9, r2 states 0.8-1.2 (about the contour). Bridge absent in r11212_s1/s2 and r222_s1
(0.31-0.37).

## Internal tree nodes (2026-11-05, node.sbatch array 54951909), same recipe as r1/r2
Path: PostProcess/branch_w_<node>/postprocess.mrc (halves Branch/<node>/h1.mrc, h2.mrc); box 360, 0.954 A/px, C2, FSC-weighted, B = 0,
unmasked, same frame. All 30 multi-leaf nodes exist: r = PostProcess/cons_w (80,598, 3.15 A), r1 (50,606, 3.18), r2 (29,992, 3.40), and:

| node | ptcls | res (A) | node | ptcls | res (A) | node | ptcls | res (A) |
|---|---|---|---|---|---|---|---|---|
| r11 | 28,823 | 3.27 | r211 | 12,170 | 3.43 | r2111 | 7,639 | 3.58 |
| r12 | 21,783 | 3.33 | r212 | 9,224 | 3.77 | r2112 | 4,531 | 3.73 |
| r21 | 21,394 | 3.37 | r221 | 5,189 | 4.09 | r2121 | 6,079 | 3.82 |
| r22 | 8,598 | 3.86 | r222 | 3,409 | 4.09 | r2122 | 3,145 | 4.29 |
| r111 | 17,911 | 3.30 | r1111 | 12,100 | 3.40 | r11211 | 4,431 | 3.82 |
| r112 | 10,912 | 3.54 | r1121 | 7,862 | 3.54 | r11212 | 3,431 | 3.95 |
| r121 | 12,189 | 3.40 | r1211 | 6,334 | 3.73 | r12111 | 3,181 | 3.95 |
| r122 | 9,594 | 3.47 | r1221 | 5,978 | 3.69 | r12112 | 3,153 | 3.95 |
| | | | r1222 | 3,616 | 3.99 | r12211 | 3,028 | 3.99 |
| | | | | | | r21111 | 4,101 | 3.82 |

Branch walk-down asked for: r -> r1 -> r11 -> r112 -> r1121 -> r11211 / r11212 (plus r12, r111).
Single-leaf nodes are the leaf map itself (PostProcess/cpose_w_<leaf>/postprocess.mrc): r1112 = r1112_g1_2, r1122 = r1122_g1_2,
r1212 = r1212_g1_2, r11111 = r11111_g1_2, r11112 = r11112_g1_2, r21112 = r21112_g1_2, r21211 = r21211_g1_2.
Leaves sit below their node (e.g. r11211_s1, r11211_s2 under r11211): PostProcess/cpose_w_<leaf>/postprocess.mrc.
