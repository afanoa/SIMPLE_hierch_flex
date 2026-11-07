---
name: klem-fortran-port-state-2026-10-16
description: klem Fortran port in progress (user: "go ahead and port so you can test faster"): worktree ~/SIMPLE_klem on batch2 (batch2 master f9bbf72ce, no state_weight_set there), new commander klem_estep registered in simple_exec_refine3D/simple_ui_refine3D, Python dump + compare harness on the toy; conventions and traps
metadata:
  type: project
---

2026-10-16 02:00-02:40. Port tree: `~/SIMPLE_klem` (git worktree, detached at batch2 master f9bbf72ce; build dir
`~/SIMPLE_klem/build`, Release/OpenMP, baseline built in 20 min on the login node; incremental = `make -j 24 install`,
a NEW source file needs `cmake .` first). Running jobs use the snapshot `~/SIMPLE_kernel/snap/5e1e6aeee`, untouched.
New code (local copy in toy55/port/): `src/main/commanders/simple/simple_commanders_klem.f90` = `commander_klem_estep`
(prg=klem_estep: projfile, dir = cell_NN_even/odd.mrc on the coarse grid, nstates, lp, infile = per-shell "shell sig2
scale", outfile = stream file with ll(K) per particle + samples; env KLEM_PINDS (rows + eo), KLEM_SS/FS/CONV/SIGN/
PHORIG/NSAMP). Registered by toy55/port/register_klem.py (exec case + UI spec via add_ui_program). First build failed
on name clashes (Fortran is case-insensitive: K vs k; `fs` shadows a core symbol) -> renamed (ncl, kk, hh, lz, ic,
fscale, sshift); rebuild running.
Reference harness: klem.py KLEM_DUMP=1 writes klem_estep_ref.npz (ll, pinds, eo, band kx/ky/wq/shell, sig2, SC, ss,
fs, conv, sign, smpd, F200/C200/R200) + cells_half/ (even/odd cell maps); toy dump = lat/klem_dump20k (20k particles,
lp 10, kmax 23, n3 47, K 51; ss -1, fs 1, conv 0, sign +1, smpd 1.464). Test = run_klem_estep_test.sbatch (waits for
make1.done, runs PHORIG 1 and 0, compare_estep.py prints corr/rel diff of V, F, C, ll and the argmax agreement).
Conventions mirrored from the prototype: half-plane band (ky outer, kx inner; kx>0 or kx==0,ky>=0; wq 2/1), F =
fftshift-centred FT x sign x exp(+2pi i (h ss fs x + k ss fs y)/box), CTF = sin(pi wl k2 (df - 0.5 wl^2 k2 cs) +
atan(A/sqrt(1-A^2))) in px units, R = Rz(-e3) Ry(e2) Rz(-e1), loc = (kx,ky,0) M (M = R or R^T), trilinear with
clamping on the centred coarse grid (phase correction exp(-2pi i (h+k+l) c3/n3) after the FFTW-style fft), ll =
-sum wq |F - C SC P|^2 / sig2 against the OTHER half's cells (eo 0 -> odd cells).
Toy trap found on the way: TOY_SMPD must be exported (view_tangent_toy defaults to 2.44 A/px) - the toy job scripts do
export it; the first dump run did not (rerun fixed).
Next: validate ll vs the dump (goal: corr > 0.999, argmax agreement ~1), then the M-step via reconstruct3D with
box_crop (or a weighted klem_mstep), the lattice smoothing / balance / readouts in a driver, then the 692k run.
Related: [[klem-kernel-latent-em-2026-10-15]], [[klem-fortran-home-state-weight-set]], [[rnasep-klem-real-data-setup]].

02:50 update: port built on the SNAPSHOT commit (worktree ~/SIMPLE_klem2 at 5e1e6aeee, build ~/SIMPLE_klem2/build;
batch2 master f9bbf72ce cannot read the newer project format). Findings from the comparison: SIMPLE's forward FFT
divides by N (2D: box^2, 3D: box^3) and already centres the phase origin (KLEM_PHORIG=0 is right; odd boxes are
unsupported: shift_phorig "even dimensions assumed", so the dump writes the coarse grid zero-padded to the full even
box, no amplitude factor). With the right phase origin but before the scale fix: V corr 1.00000, F corr 1.00000, C
exact, ll per-particle centred corr 0.95, argmax agreement 0.41 -> scales applied (F x box^2, V x nc^3), rebuilt
(make3.done), test rerunning. Dump = lat/klem_dump20k (job 53904140).

02:55 E-STEP VALIDATED: prg=klem_estep (build ~/SIMPLE_klem2/build) vs the Python dump on 20k toy particles x 51
cells at lp 10: V rel diff 0.0000, F 0.0000, C exact, ll per-particle corr 1.0000 (rel diff 0.0000), argmax
agreement 1.000, responsibilities identical; with KLEM_PHORIG=0, F x box^2, V x nc^3. ~40 s for 20k x 51 on 16
threads (409k ~14 min single node; ~1.5 min over 10 parts). Next: M-step (the Python M-step already runs a Fortran
kernel via ctypes: port that source), then the board driver (lattice smoothing, balance, bands, readouts, junk).

03:20 BOARD DRIVER WRITTEN: src/main/commanders/simple/simple_klem_core.f90 (module simple_klem_core: band, particle
loading, projection, cross-half E-step, nearest-voxel M-step, lattice kernel, Sinkhorn, full-box volume writer) and
simple_commanders_klem_board.f90 (prg=klem_board: projfile, infile = config "key value" lines (n1 n2 nfree lps sigmas
shrink eps0 balance balance_until focus_mask subset nsub seed ss fs conv sign ctfabs write_halves), dir = output;
random init, per-band consensus + per-shell scale/noise calibration, kernel-smoothed Wiener M-step with shrinkage,
focus via the full even box, softmax + Sinkhorn balance, readouts; outputs klem_latent.txt, klem_states.txt,
klem_resp.bin, klem_halves.bin, cells/*.mrc at the full box). klem_fortran2npz.py converts the .bin files to the npz
layout the Python tools read. Registered via register_klem.py (klem_estep + klem_board). Build make4 running; test =
run_klem_board_toy.sbatch (20k toy particles, lps 12 12 10, no balance, same conventions as the Python dump board;
readouts: state_coverage.py vs truth for both boards, deposited_match.py on the Fortran cells). Not ported (v1): pose
marginalisation, averages/unmix modes, pca/npz init, per-particle weights (pw = 1).

03:50 BOARD VALIDATED ON THE TOY: prg=klem_board (20k particles, lps 12 12 10, random init, no balance) = 99 s on 16
threads; truth purity 0.17 (5/10 own a cell) vs the Python reference 0.16 (4/10) on the same particles/schedule;
free mass 0.18 (toy junk 20 %). Per-round cost ~3 s for E+M at 20k x 51; ~43 s per band for loading + calibration
-> next: load once at the finest band and mask coefficients per band (3x fewer loads). The Python tooling
(klem_fortran2npz.py -> klem_resp.npz / klem_halves.npz; cells/*.mrc) consumes the outputs unchanged, so
expand_families.py / node_states.py / derive_regions.py / rank_by_target.py / deposited_match.py run on Fortran boards.
04:00: load-once-per-schedule added (klem_estep_ll / klem_mstep_accum / klem_project_cell take kc = the current band
limit and skip coefficients beyond it; the board loads at min(lps)). First real-data run: closed_focal_fortran_rnasep.sbatch
(family 1, 409k, focal mask, config written from the Python run's calibration line, outputs in
klem_fb/closed_focal_fortran/board + the same cell-scoring/reconstruction tail as closed_focal_rnasep.sbatch).
04:40: the whole-reference local pose step ported (klem_pose_step in simple_klem_core.f90: 27-rotation grid at step and
step/2, Gaussian prior, posterior-mean update, against the argmax cell of the other half; klem_rotmat2euler back to
SIMPLE angles; board config keys pose_every/pose_step/pose_prior; output klem_poses.txt "pind e1 e2 e3"). Written in
toy55/port/, NOT yet built: do not `make install` in ~/SIMPLE_klem2 while a job runs from that build (the real-data
focal board job). Next: build after it ends, test on the wide toy (truth: lever 3.99 -> 2.89 deg expected), then
the real docked 10k (expected contact 1.85, FSC 6.9 A).
04:55 FIRST REAL-DATA RUN (53904203, family 1, 409k, focal mask, 9 rounds, random init, balance 0.5): 41 min on 16
threads (loading 1008 s; rounds 72-240 s; cells written) vs the Python's 2 h. But weaker: docked corner max +0.59
(Python +1.11), focal map 1 density 89 % (102 %), cross-half reproducibility median -0.28 (0.95): noise-driven cells at
the fine bands -> port devw (deviation-FSC shell weights, sig2 = sig2_base / clip(2F/(1+F), devw, 1)), pnoise
(per-particle noise factor) and the npz/latent init before judging. Loading: parallelise the reads next.
05:05: board v2 additions written and building (make10): devw (deviation-FSC shell weights from round 2), pnoise
(per-particle noise factor from the residual against the other half's consensus), init_file (lattice positions from
klem_init_positions.py: standardised latent PCs 1-2, 5-95 pct onto the grid; the meta's z692), flush(logfhandle).
Real-data v2 run = closed_focal_fortran_rnasep.sbatch ZTAG=_v2 (DEVW 0.05, PNOISE 1, INITFILE) to compare with the
Python board (docked corner +1.11, reproducibility 0.95) and with v1 (+0.59, -0.28). Toy pose-step test
(run_klem_board_pose_toyw.sbatch, projw.simple) running first. Build make9 (pose step) = OK.
05:15 POSE STEP VALIDATED on the wide toy (53904329): Fortran board, true conformer-18 particles, pose_every 2 /
step 3 / prior 4: pose error 3.81 -> 2.90 deg (p90 7.01 -> 5.68), arm 0.80 -> 0.82; Python gave 3.99 -> 2.89. 32 s
for 10k particles x 6 rounds incl. 3 pose steps on 16 threads.
06:40 V2 REAL-DATA RESULT (53904364: devw 0.05, pnoise, latent init, balance 0.5, focal mask, 409k, 9 rounds): 31 min
on 16 threads (Python 2 h). Cross-half reproducibility median 0.99 (49/49; v1 -0.28; Python 0.95); the docked corner
is the SAME cells as the Python board (01_04, 00_04, 00_05, 01_05, ...), scores +0.64 (Python +1.11: the Fortran cells
carry smaller deviations, no pose marginalisation); the top-3-cell docked set (26.7k) reconstructs identically to the
Python's (corr +0.42 / -0.50, density 98 % vs +0.41 / -0.51 / 102 %). The Fortran board is a drop-in for the real-data
recipe at 4x the speed; next: parallel particle reads (17 min of the 31), then the tree driver around it.

11:15 LOADER: klem_load_particles now sorts the subset into stack order (map_ptcl_ind2stk_ind + hpsort on stkind*1e7+ind)
and reads through the buffered stack_io reader (one open per stack, sequential) instead of image%read per particle
(open/close + random seek per particle). Validated on the 20k toy: argmax identical for all 20k, z within 7e-4,
56 -> 47 s on a single sequential stack; the gain on the real 10-stack permuted project is to be measured. Built in
~/SIMPLE_klem2/build (make11); the running f692 chain and the toy board jobs use frozen copies (rnasep_klem/f692/klem2_bin).

13:30 DISTRIBUTED BOARD VALIDATED: prg=klem_board nparts=N (master, no particles) + prg=klem_board_part (persistent
SLURM workers, one job each, file exchange in DIR/distr: CONS/SCALE/NOISE/MSTEP/ESTEP/FINISH requests, half-sphere
packed cubes, per-part Sinkhorn balance), source port/simple_klem_distr.f90 + dispatch in the board commander, UI
'nparts' on klem_board, worker program registered by register_klem.py. Toy 20k: cells corr 0.99 with the sequential
board (two different cells 0.95), truth readout identical, argmax agreement 0.71 (calibration subsample differs).
Real family-1 feedback board 307k: 442 s on 6 workers x 8 threads vs 1012 s sequential at 16 threads, reproducibility
0.98 (49/49). Launch traps fixed on the way: qsys submits "<cwd>/<script>" (use relative script names), the default
worker executable is simple_private_exec (pass exec_bin='simple_exec'), a multi-command script runs its jobs
sequentially (one exec_simple_prg_in_queue_async per worker), workers read DIR/distr, band-change reallocation.
Wrapper: fboard.sh NPARTS=N (+ NTHR per worker); next-run chain f692_chain5.sbatch (NPARTS 10, NTHR 8, concurrent
families, parallel regions, NP refine3D_states parts). WAIT_MAX 1 h frees workers if the master dies.
