---
name: flex2d-branch-and-test-setup
description: flex2D lives on branch flex2D (Mac + batch2 worktree ~/SIMPLE_flex2D); how it is built, run and scored on the 10028 8000-particle subset
metadata:
  type: project
---
Branch `flex2D` (from master d5e61b7fc, 2026-09-21). Mac repo pushes to the batch2 remote through the
scratch ref `flex2D-in` (`git push batch2 flex2D:refs/heads/flex2D-in`, the branch itself is checked out
in the worktree there); on batch2: `cd ~/SIMPLE_flex2D && git reset --hard flex2D-in && bash ~/build_flex2D.sh`
(incremental, ~1 min for flex2D-only edits, ~3 min when the pftc parent module changes; done-file
`build/make.done`, binary `~/SIMPLE_flex2D/build/bin/simple_exec`).
Tests: `/mnt/beegfs/elmlund/afan/10028/flex2D_tests/` — `ref.simple` (the reference 2D run, 300 classes,
105k ptcls, box 360, smpd 1.34), `ref_ptcl2D.txt` (its ptcl2D dump), `compare_poses.py`, `clsres.py`,
`mrc_montage.py`; per-test dirs tN/ with `run_srun.sh` (selection nran=8000 -> flex2D -> dumps -> compare).
Run interactively on a compute node: `srun -p csbdevel -c 48 --mem=180G -t 60 bash tN/run_srun.sh`
(numeric SIMPLE calls SIGILL on the login node; print_project_field survives there). One test = ~90 s.
Scoring: within-reference-class circular spread of (e3_test - e3_ref) (random = 104 deg), shift residual
rms, NMI/purity of labels, per-class FRC resolution from cls2D (`clsres.py`). Reference restored from an
8000 subset at box 182 = 10.5 A median (the fair baseline; the full run is 5.2 A).
Doc: `doc/implementation_notes/flex2D_proposal_2026_09_21.md`. See [[flex2d-findings-2026-09-21]].
