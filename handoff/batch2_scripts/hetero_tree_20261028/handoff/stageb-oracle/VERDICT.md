# What limits stage B's state assignment: poses, references or noise

Oracle refine=state_eval arms, 2026-11-06/07. Side session "stage-B oracle" of the hetero_tree paper thread.

## Verdict

**The references limit stage B: their number at the root and their purity inside a split. The poses do not. Noise sets the
ceiling, and stage B reaches it only where its references already match the classes.**

1. **Poses are not the limit.**
   - Stage-B poses are within 2.2-2.6 deg (median, non-junk) of the truth on qtoy2 and its children. Fewer than 1 % are off by
     more than 10 deg.
   - True poses with stage B's own maps change no accuracy: qtoy2 root group 0.951 -> 0.947; r1 class 0.756 -> 0.758; r4
     0.912 -> 0.904.
   - True poses do move labels (root 11 %, r1 18 %), but they fix as many as they break (root: 22 % of the misassigned
     fixed, 12 % of the right ones broken). Stage B's poses and maps are co-adapted.
   - Poses pay only once the references are pure: +1.5 to +4 points (root OB -> OT +0.017; r1 OB -> OT +0.041, CB -> CT
     +0.025).
   - The exception is qtoyflex, with a heavy pose tail (median 4.3 deg, 10.7 % above 10 deg, 4.5 % above 30 deg). There, true
     poses alone lift holo F1 0.617 -> 0.682.

2. **The references are the limit, node by node.**
   - **Root, too few references.** qtoy2's root has K = 5 for 8 classes. Its 5 references are pure as groups: 95.1 % of the
     non-junk particles sit in their group's state. Pure references at K = 5 give 95.3 %; with true poses as well (the K = 5
     ceiling), 96.3 %. So at its own K the root is 1.2 points from the ceiling. But K = 5 merges 2|1, 18|22 and 43|49. One
     reference per true class at the same stage-B poses lifts class accuracy 0.674 -> 0.919 (ceiling 0.936) and closed F1
     0.55 -> 0.95.
   - **r1 (classes 1 vs 2, K = 2), impure references.** Stage B's two references are 80/20 mixtures of the two classes. Pure
     references at the same K and poses: 0.756 -> 0.839. With true poses: 0.864. Ceiling 0.891. Pure references close 62 %
     of the gap at stage-B poses. The 80/20 pair is a fixed point: state_eval at stage B's maps reproduces it.
   - **r4 (closed 18 vs 22, K = 2), at the ceiling.** 98.9 % of closed and 97.9 % of class 22 are in their own state.
     Pure references and true poses change nothing. Closed F1 (0.896) is held down by junk, 17 % of the node. A junk
     reference lifts it to 0.939 at the same poses.
   - **qtoyflex root, no closed reference.** K = 3 is spent on two arm halves plus holo, so closed F1 is 0.168. With
     references for arm, closed and holo at the same K and the same stage-B poses: 0.455 (recall 0.88). With true poses:
     0.570; ceiling 0.634. Closed precision is capped by arm conformers next to closed: 8 % of the arm goes to the closed
     reference even at true poses.

3. **Noise.**
   - The ceiling arm (true poses, references rebuilt from the truth at true poses) leaves 6.4 % class error at the qtoy2 root
     (K = 9), 10.9 % in r1, and closed F1 0.63 on qtoyflex.
   - The ceiling's errors fall with snrq: r1 0.84 in the lowest quintile vs 0.93 in the highest.
   - Stage B's extra errors do not. The r1 gap to the ceiling is 12-15 points in every snrq quintile. Where stage B sits at
     its ceiling (root at K = 5, r4), its leftover errors are the low-snrq ones.

4. **The "32 % confidently wrong" at the qtoy2 root** (mu_qtoy2's stage-B posteriors; mm2 by hard labels gives the same):
   - 32.5 % wrong = 27.8 % **merged** + 4.7 % misassigned. Merged means classes 2, 18 and 43, whose majority state belongs to 1,
     22 and 49. Their class owns no state at K = 5, which is why their posterior mass on the right class is exactly 0.
   - Of the 13,294 confidently wrong (max posterior >= .99), 96.7 % are merged. Against the confidently right they are
     indistinguishable: pose error median 2.30 vs 2.49 deg, snrq median 1.05 vs 1.03, B-factor 98 vs 99, no junk. Junk itself
     is mostly unconfident (mean max posterior 0.78; 12.5 % >= .99).
   - One reference per class at the same poses puts 92.7 % of the merged in their own class.
   - The truly misassigned 4.7 % are another population:
     - mostly unconfident (max posterior median 0.87; 17 % >= .99; mean mass 0.11 on their own class's state);
     - lower quality (snrq median 0.86 vs 1.01; B-factor 107 vs 99);
     - a pose tail (median 3.1 deg; 2.8 % above 30 deg; the misassigned share rises from 2 % at < 1 deg to 10 % at 5-10 deg).
     - AUROC against the right ones: pose error 0.62, snrq 0.61.
     - Pure references fix 41 % of them at stage-B poses and 82 % at true poses.

5. **Where gains remain: in the references.**
   - The tree's final qtoy2 leaves (best-leaf F1 1: 0.77, 2: 0.69, 18: 0.89, 22: 0.73, 28: 0.94, 43: 0.85, 49: 0.83,
     101: 0.94) are below one state_eval pass with one pure reference per class at the root's stage-B poses (0.88, 0.84,
     0.95, 0.87, 0.96, 0.90, 0.90, 0.97). The largest gaps are the r1 pair (1, 2) and 22.
   - What to improve:
     - the split that should separate a close pair (r1's 80/20 fixed point);
     - a junk reference wherever junk sits next to a small state (r4 closed; the 22 leaf is 19 % junk);
     - for a continuum, a root state for the minority conformer rather than a second arm state (qtoyflex).
   - Better poses are worth 1.5-4 points after the references are fixed, and more under heavy pose tails (qtoyflex holo).
     They are worth nothing before.
   - This agrees with the day's earlier verdicts: stage B is seed-invariant, and the objective is sufficient once references
     are good, so the work is in getting the references right (the split), not in the posterior or the poses.

## Numbers

Accuracies are over non-junk particles. **class** = in the state of their own class: the reference's class for oracle arms,
otherwise the state's majority class. **group** = in the state of their coarse group (stage B's grouping). Poses: B = stage B,
T = truth carried into the stage-B frame. Maps: SB = stage B's; Ctl = rebuilt from stage B's labels; coarse = the truth merged
into stage B's groups; oracle = one per truth class (junk is a class). The suffix B/T says at which poses a map set was rebuilt.

qtoy2 root (node_r; 60,000 = 54,000 + 6,000 junk; stage B K 5; pose error 2.63 deg):

| arm | poses x maps | K | NMI all | class | group | F1 closed 18 | F1 holo 101 | F1 junk |
|---|---|---|---|---|---|---|---|---|
| P | B x SB | 5 | 0.662 | 0.674 | 0.951 | 0.550 | 0.943 | 0.228 |
| C | B x Ctl | 5 | 0.662 | 0.674 | 0.951 | 0.550 | 0.943 | 0.228 |
| CB | B x coarse-B | 5 | 0.667 | 0.677 | 0.953 | 0.553 | 0.945 | 0.275 |
| TB | T x SB | 5 | 0.660 | 0.670 | 0.947 | 0.572 | 0.943 | 0.243 |
| CT | T x coarse-T | 5 | 0.699 | 0.684 | 0.963 | 0.553 | 0.962 | 0.386 |
| OB | B x oracle-B | 9 | 0.786 | 0.919 | 0.965 | 0.947 | 0.965 | 0.776 |
| OBT | B x oracle-T | 9 | 0.749 | 0.891 | 0.954 | 0.928 | 0.963 | 0.744 |
| OTB | T x oracle-B | 9 | 0.803 | 0.903 | 0.961 | 0.963 | 0.983 | 0.896 |
| OT | T x oracle-T (ceiling) | 9 | 0.856 | 0.936 | 0.976 | 0.965 | 0.993 | 0.932 |

The NMI ceiling of a perfect assignment into stage B's 5 groups is 0.845. Per-class F1, P / OB / OT: 1 .69/.88/.90,
2 .52/.84/.87, 22 .52/.87/.92, 28 .94/.96/.98, 43 .46/.90/.94, 49 .65/.90/.94.

qtoy2 r1 (classes 1 vs 2, plus 22 and junk; 21,951; K 2; pose error 2.19 deg; frame rotated 1.12 deg from the truth):

| arm | K | class | group | F1 1 | F1 2 |
|---|---|---|---|---|---|
| P | 2 | 0.756 | 0.796 | 0.789 | 0.720 |
| C | 2 | 0.756 | 0.797 | 0.790 | 0.720 |
| CB | 2 | 0.839 | 0.878 | 0.862 | 0.815 |
| TB | 2 | 0.758 | 0.794 | 0.795 | 0.717 |
| CT | 2 | 0.864 | 0.904 | 0.888 | 0.840 |
| OB | 4 | 0.850 | 0.887 | 0.884 | 0.857 |
| OT | 4 | 0.891 | 0.913 | 0.917 | 0.894 |

Class accuracy per snrq quintile, P / OT: .70/.84, .74/.87, .77/.89, .78/.92, .79/.93.

qtoy2 r4 (closed 18 vs 22, plus strays and 17 % junk; 9,413; K 2; pose error 2.23 deg):

| arm | class | group | F1 closed | F1 22 | F1 junk |
|---|---|---|---|---|---|
| P | 0.912 | 0.974 | 0.896 | 0.823 | 0.286 |
| CB | 0.913 | 0.976 | 0.905 | 0.818 | 0.301 |
| TB | 0.904 | 0.970 | 0.910 | 0.798 | 0.324 |
| CT | 0.905 | 0.972 | 0.924 | 0.791 | 0.346 |
| OB, K 4 | 0.869 | 0.972 | 0.939 | 0.875 | 0.708 |
| OT, K 4 | 0.864 | 0.969 | 0.955 | 0.881 | 0.818 |

OB and OT class accuracy falls because the 328-particle class-2 reference draws 10 % of class 22.

qtoyflex root (kind truth: arm continuum 44,366, closed 3,600, holo 6,016, junk 6,000; stage B K 3 = two arm states + holo;
pose error 4.34 deg, 10.7 % > 10 deg):

| arm | K | NMI all | class | group | F1 closed | F1 holo | F1 junk |
|---|---|---|---|---|---|---|---|
| P | 3 | 0.167 | 0.848 | 0.905 | 0.168 | 0.617 | 0.176 |
| C | 3 | 0.167 | 0.847 | 0.904 | 0.168 | 0.616 | 0.177 |
| TB | 3 | 0.214 | 0.874 | 0.929 | 0.166 | 0.682 | 0.184 |
| NB (arm/closed/holo, B) | 3 | 0.269 | 0.806 | 0.916 | 0.455 | 0.650 | 0.163 |
| NT (arm/closed/holo, T) | 3 | 0.405 | 0.888 | 0.959 | 0.570 | 0.780 | 0.171 |
| OB (+ junk reference) | 4 | 0.296 | 0.696 | 0.923 | 0.501 | 0.685 | 0.455 |
| OT (ceiling) | 4 | 0.538 | 0.850 | 0.961 | 0.634 | 0.838 | 0.764 |

The coarse arms collapse to K 2 here (stage B has no closed state), so NB/NT are the K-matched oracle.

Wrong particles, qtoy2 root (`wrong_mu.txt`; mm2 in `wrong.txt` agrees within 0.002):

| group | n | max posterior >= .99 | pose error median / >30 deg | shift median | snrq median | B-factor |
|---|---|---|---|---|---|---|
| right | 36,473 | 0.88 | 2.50 / 0.0 % | 0.33 px | 1.006 | 99 |
| merged (K = 5) | 14,985 | 0.86 | 2.33 / 0.0 % | 0.32 px | 1.013 | 100 |
| misassigned | 2,542 | 0.17 | 3.06 / 2.8 % | 0.40 px | 0.862 | 107 |
| junk | 6,000 | 0.13 | 10.8 / 23.9 % | 2.10 px | 0.995 | 100 |

## Method

- **Arms.** Every arm is the node's own stage-B refine3D call (`stageB/hetero_step.sh`) at `refine=state_eval maxits=1`,
  with the state_eval build and `SIMPLE_SEED=1`, on a copy of `stageB/b.simple`.
  - Only the ptcl3D slots e1 e2 e3 x y and state change. `orc_patch.py` byte-patches the binoris records and checks them
    against print_project_field before and after.
  - state_eval recomputes the grid direction, the in-plane index and the shift from those slots.
  - The noise weights are global per shell and half, with no per-state term, so starting labels give no arm an advantage.
- **True poses.** From the simulation's `all.txt` (x = y = 0), carried into the stage-B frame by a fitted global rotation and
  3D shift (`orc_prep.py`).
  - Roots: rotation 0.02-0.21 deg, shift 1.87 px (rms shift error 1.68 -> 0.43 px).
  - Children: 1.1-1.3 deg rotation.
  - The direction stays on the nspace 5000 grid in every arm, true poses included.
- **Oracle references.** Rebuilt by the same refine3D call (PCG + pcg_focus + ML regularisation) from fixed labels.
  - `SIMPLE_STATE_PRIOR` holds 0 on the wanted label and -1e4 elsewhere, at maxits 2. The first references come from the
    node's recB reconstruct3D call. Every forced run kept 100 % of its labels.
  - The control C (stage B's own labels and poses) matches production to 0.001, which validates this route against stage
    B's maps.
- **Label sets.**
  - Truth classes (junk one class; in a child node, classes under max(300, 2 %) join junk).
  - Coarse: each class merged into its majority stage-B state, K <= K_B.
  - No-junk (`orc_nj.sbatch`): junk is left out of the maps and assigned like everyone else.
- **Wrong particles** follow mest_post.py: a state's owner is its majority non-junk class. "Merged" means a particle sits in
  its own class's majority state, but another class owns that state. Posteriors come from mu_qtoy2's
  `stageB/prob_state_evidence.txt`. The mm2 build writes no evidence, so mm2 is scored by hard labels.

## Where things are

- **Scripts** (copies next to this file; the working copies are in `math_reworks/latent/`):
  - `orc_prep.py`: labels, priors, true poses in the stage-B frame, pose error.
  - `orc_patch.py`, `orc_verify.py`: the slot patch and its check.
  - `orc_score.py`, `orc_wrong.py`: readouts.
  - `orc_run.sbatch` (first version), `orc_run2.sbatch` (resumable, adds CT), `orc_nj.sbatch` (the no-junk oracle).
- **Runs**: `math_reworks/latent/oracle/<key>/`, with key qtoy2 (its `mu/` holds the mu_qtoy2 root prep), qflex, qtoy2_r1,
  qtoy2_r4. Each has:
  - `prep/` (`prep_report.txt` = frame fit and pose error);
  - `R_*` (recB first references), `F_*` (forced-label runs whose `recvol_stateNN.mrc` are the oracle maps);
  - `A_*` (arms; `after.txt` = the arm's ptcl3D dump);
  - `score.txt`, `wrong.txt`, and `score_nj.txt` (qflex only), `wrong_mu.txt` (qtoy2 only).
- **Jobs** (logs `math_reworks/latent/logs/orc*_<id>.out`): 55010527 qtoy2, 55010528 qflex, 55010712 r1, 55011776 r4,
  55012117 qflex no-junk, 55012419 qtoy2 + CT. All finished with exit 0.
- **Resume or extend:** `sbatch -J orc_<key> --export=ALL,KEY=<key>,NODE=<gate node>,TRUTH=<truth file>,ALLTXT=<all.txt>[,ORC_MINN=300,ORC_MINFRAC=0.02] orc_run2.sbatch`.
  It skips finished steps. A child node needs ORC_MINN and ORC_MINFRAC. About 40 CPUs at peak and 20-30 min per node.
- **Follow-up worth a run:** the same arms on the other children (r2 = 43|49, r3, r5) and on r131's closed nodes. And r1's
  split from pure references: does stage B, started at CB's maps, keep them, or fall back to 80/20?
