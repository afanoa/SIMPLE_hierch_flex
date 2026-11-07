---
name: sieve-sigma2-merge-bugs
description: particle_sieving dies at its merges from two sigma2 bugs; the intent(out) aliasing one is fixed by a local patch, and finished chunks can be re-merged without redoing compute
metadata:
  type: project
---

Hit both on 2026-09-18 running `particle_sieving` on trpm4 and flipqr.

**Bug 1 (fixed upstream).** Coarse->fine promotion called `merge_chunk_projfiles`,
which needs each chunk's `sigma2_state.bin` — but `cleanup_chunk` deletes it.
Crash: `sigma2 state file does not exist ... simple_projfile_utils.f90:337`.
Fixed by `bc63531a0` (2026-09-16), which switches that call to
`merge_chunk_projfiles_without_sigma2`. A build older than that WILL die ~5 min in.

**Bug 2 (still in master at 34e69f0c4; patched locally).** The final combine dies with
`invalid sigma2 state grouping policy; simple_projfile_utils.f90:306` **after all
chunk work finishes**. Root cause is an aliasing violation:
`sigma2_state_init_header` declares `type(sigma2_state_header), intent(out) :: header`,
and the caller passes `target_header` as that argument *while also passing its own
components* (`target_header%kfrom`, `%box`, `%smpd`, `%grouping`) as inputs. The
`intent(out)` makes the variable undefined on entry, so `grouping` arrives as 0 and
fails the `GLOBAL(1)|STACK(2)` check. **Two sites: lines 302 and 850.** Fix = copy the
fields into locals before the call. Patch both; 850 is a different merge path.

**Recovering a run that died at the final merge.** The sieve has NO resume logic, so
re-running repeats hours of chunk compute. Instead the finished chunks can be
re-merged directly: the completed fine-chunk projects live at
`chunks_fine/chunk_fine_<N>/chunk_fine_<N>.simple` (with their `sigma2_state.bin`
intact — only *coarse* chunks get theirs cleaned), and `combine_completed_chunks`
is just `merge_chunk_projfiles(projfiles, completedir, proj, write_proj=.false.)`
followed by `proj%write()`. A ~40-line program doing that recovered both datasets in
minutes. Verify by comparing `state=1` count against the log's
`TOTAL NUMBER PARTICLES ACCEPTED` — they matched exactly.

**Compiling a one-off against the SIMPLE build** (both tools needed this):
```
gfortran -O1 -fopenmp -o tool tool.f90 -I$B/modules -I$B/lib/simple \
  $B/lib64/libSIMPLE3.0.0.a /usr/lib64/libarpack.so /usr/lib64/libopenblas.so \
  -lstdc++ -lz -ljpeg -ltiff -lfftw3f -lfftw3f_threads -lfftw3 -lm -ldl -lpthread
```
arpack + openblas are both required (`ssaupd_`, `sgemm_`). Project paths are
`class(string)`, not `character` — wrap with `string(trim(x))`.

Counting selected particles: grep `" state="` WITH the leading space, or `mi_state=0`
inflates the state=0 count by one per particle.

See [[particle-sieving-job-model]], [[trpm4-flipqr-sieve-campaign]].
