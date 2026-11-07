---
name: klem-f692e-branch-snapshot
description: "~/SIMPLE_f692e: tag f692e-2026-10-18 = exact f692e code; branch klem-consolidated = prg=klem_tree, the whole route in Fortran, pushed to the PUBLIC fork afanoa/SIMPLE_hierch_flex; its full run f692k reproduces f692r at result level (2026-10-20)"
metadata:
  node_type: memory
  type: reference
  originSessionId: aa724bb9-dd55-4c19-9967-ca5632aea4b5
  modified: 2026-10-20T23:56:25.634Z
---

Made 2026-10-19 at the user's request ("preserve the code from that Oct 18 run at 2:05 pm ... retain the current code too").
Before this NOTHING of klem was under version control: the Fortran lived untracked in batch2 ~/SIMPLE_klem2 (detached at
5e1e6aeee + 3 modified tracked files) and the pipeline in untracked toy55/ (local) and routes/toy55 (batch2, 66 extra files).
- Tag f692e-2026-10-18 (5ca3c3b7b): src/main/commanders/simple/simple_klem_{core,distr,tree}.f90 + simple_commanders_klem{,_board,_tree}.f90
  as built into rnasep_klem/f692e/klem2_bin (2026-10-17 18:36); toy55 code as of 2026-10-18 14:05 (306 unchanged + 5
  restored: f692_family.sh, consolidate.sbatch, consolidate_states.py, merge_data_derived.py, pooled_dd.sbatch); README
  toy55/F692E_RUN.md (launch command, downstream steps, outputs, method).
- Head 59a1e1b28: current klem sources (batch2, with gridding/specmatch/hp off by default) + all current scripts.
**How to apply:** to rerun f692e exactly, build from the tag (register the klem commanders as in ~/SIMPLE_klem2) or
use the frozen binaries in rnasep_klem/f692e/klem2_bin; never push without the user's permission. Method for future
reconstructions: undo recorded python/sed edits from the session jsonl and replay them to verify (scratchpad recon/undo_patches.py).
VALIDATED 2026-10-20 10:00 (f692e_repro_test.sbatch, toy55): the Oct 17 smoke parent board (f692smoke7/parent: 100k
real RNase P particles, 7x7 + 2 free, 20/14/10 A, 6 local workers x 8 threads; made with the SAME binaries the f692e run
froze, md5 58154f54e86c) rerun on (a) the frozen f692e klem2_bin and (b) the current ~/SIMPLE_klem2/build (Oct 19 12:50,
md5 ccf5d3330eb2): both BIT-IDENTICAL to Oct 17 (100,000 latent rows, 51 cell maps, max diff 0) -> the board is
deterministic and the current build with default settings (gridding/specmatch/hp off) = the f692e board. This is the
small regression test: rerun it after any klem change (~3 min per arm; private project copy + update_project qsys local).
Not covered: the f692e pipeline scripts downstream of the board (tag versions differ from current) - a full f692e
rerun from the tag would take ~20 h.
FORK (user 2026-10-20): https://github.com/afanoa/SIMPLE_hierch_flex = remote `fork` (git@github.com:afanoa/SIMPLE_hierch_flex.git,
read access OK; its master abeeebfba). Plan agreed: once run 1 (f692r) confirms the Oct 18 result, consolidate on klem-f692e
(klem commanders + pipeline scripts with the 2026-10-20 fixes + regression test + monitor; keep the tag) and push THAT
BRANCH (+ tag) to `fork` - never to origin (hael/SIMPLE); show the commit summary and get the go before pushing.
REPRODUCED 2026-10-20 14:10 (user: "f692r looks good, contains the closed state"): f692r = the frozen code rerun
(Oct 20 10:22-13:00, rnasep_klem/f692r; batch2 routes/toy55_f692e + the 2026-10-20 fixes). Consolidation of the build
(less bash/python glue) + push of branch and tag to `fork` handed to a SEPARATE session at the user's request.
CONSOLIDATED 2026-10-20 (session happy-borg): branch klem-consolidated (on klem-f692e) = prg=klem_tree via simple_private_exec;
no bash/Python left (doc/klem/f692e_route.md, launcher scripts/klem/klem_tree.sbatch, config scripts/klem/rnasep_f692.txt,
board regression scripts/klem/klem_board_repro.sbatch -> VERDICT IDENTICAL). Pushed to `fork` at b47045102 + the tag; the
fork is a PUBLIC fork of public hael/SIMPLE (only the PR banner is private) - user chose to leave it public.
Full run f692k (rnasep_klem/f692k, 3.5 h on norm) vs f692r: families identical; landscape 78/110/154/256 vs 68/96/149/246
states at sim 0.6-0.9 with the same 9/9/10/10 conformers identified; region pool 67 vs 68 states, strict recovery 5 vs 5
of 55; zoom pool 234 vs 218, recovery 12 vs 11; arm readout 0 closed in both. Region boards differ run to run
(refinement randomness, also f692r vs f692e), so compare runs at result level. Readouts are twins of the Python (0 of 301
states differ). Resume a run on a newer build: remove run/klem_bin/bin; the launcher refuses a run dir it did not start.
