#!/bin/bash
#SBATCH -J klem_panel
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# one benchmark set through the general pipeline, frozen settings (panel_run.sbatch with board seeds): (0) whole-particle klem board
# (7 x 7, no free component, 20 -> 10 A; reused from ff2_<TAG> when present) -> its cells in K spectral clusters = the seeds;
# (1) discovery on all particles = K-state refine3D_states from those seeds (pose_policy local, 20 -> 12 A, 6 its); (2) leaf tree from
# the leaves' half-maps (state_merge.py, complete linkage, 40-10 A); (3) merge-up along the tree (bu_merge.sbatch: splitbench9 on the
# even/odd rows, cross-half EM, view x defocus stratified null); scorecard after every stage.
# env TAG (sim<TAG>, proj<TAG>, lat_<TAG>), K (20), OUT (rnasep_klem/panel2_<TAG>); scorecard GT/CLASSES/CODES (environment)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; V=/mnt/beegfs/elmlund/afan/routes/toy55_f692v2; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
TAG=${TAG:?}; K=${K:-20}; O=${OUT:-$K_/panel2_$TAG}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') panel2 $TAG] $*"; }; GT=${GT:-$W/sim$TAG/gt_particles.txt}
NR=$(grep -c . $W/lat_$TAG/ptcl3D_r0.txt)
# ---- (0) board seeds
BD=$K_/ff2_$TAG/whole/board; [ -s $BD/klem_latent.txt ] || BD=$O/whole/board
if [ ! -s $BD/klem_latent.txt ]; then
  export FB=$K_/f692e/klem2_bin CTFABS=0 PSCALE=0 SERVER_LIST= NPARTS=; mkdir -p $O/whole; seq 1 $NR > $O/whole/all_rows.txt
  [ -s $O/meta_toy.npz ] || $PY $W/make_toy_meta.py $W/lat_$TAG/ptcl3D_r0.txt $O/meta_toy.npz > /dev/null
  cp -f $W/proj$TAG/proj$TAG.simple $O/whole/p.simple; log "whole-particle board on $NR particles"
  NTHR=64 N1=7 N2=7 NFREE=0 FOCUSMASK= META=$O/meta_toy.npz bash $V/fboard.sh $BD $O/whole/p.simple $O/whole/all_rows.txt "20 14 12 10 10" "1.5 1.2 1.0 1.0 1.0" > $O/whole/board.out 2>&1 || { log "board FAILED: $(tail -1 $O/whole/board.out)"; exit 1; }
fi
[ -s seed_states.txt ] || $PY - <<PYEOF
import glob, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23])).reshape(n, n, n).astype(float)
fs = sorted(glob.glob('$BD/cells/cell_*.mrc')); V = np.array([rd(f) for f in fs]); D = V - V.mean(0); X = np.array([np.fft.fftn(d).ravel() for d in D]); X /= np.linalg.norm(X, axis=1)[:, None]
S = np.maximum(np.real(X @ X.conj().T), 0); np.fill_diagonal(S, 0); dg = S.sum(1) + 1e-9; ev, U = np.linalg.eigh(np.eye(len(S)) - S / np.sqrt(np.outer(dg, dg))); k = min($K, len(fs))
E = U[:, :k] / np.maximum(np.linalg.norm(U[:, :k], axis=1)[:, None], 1e-12); rng = np.random.default_rng(0); best = None
for rep in range(20):
    C = E[rng.choice(len(E), k, replace=False)]
    for it in range(100): lab = np.argmin(((E[:, None] - C[None]) ** 2).sum(-1), 1); C = np.array([E[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
    sse = ((E - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
lab = best[1]; L_ = np.loadtxt('$BD/klem_latent.txt', comments='#'); pind, cell = L_[:, 0].astype(int), L_[:, 3].astype(int); ncell = len(fs)
out = np.zeros($NR, int); on = (cell >= 1) & (cell <= ncell); out[pind[on] - 1] = lab[cell[on] - 1] + 1
np.savetxt('seed_states.txt', out, fmt='%d'); print('board seeds: %d cells in %d spectral clusters, %d of %d particles' % (ncell, k, (out > 0).sum(), $NR))
PYEOF
# ---- (1) discovery
if [ ! -s refined_states.txt ]; then
  K=$(awk '$1 > 0' seed_states.txt | sort -n | uniq | wc -l); log "discovery: $K states from the board seeds"
  cp -f $W/proj$TAG/proj$TAG.simple p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local pgrp=c1 mskdiam=220 lpstart=20 lpstop=12 maxits=6 ml_reg=no nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR' ref.log | cut -c1-150)"; exit 1; }
  RD=$(ls -td $O/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > refined_states.txt
fi
RD=$(ls -td $O/*_refine3D_states | head -1); sc(){ $PY $W/scorecard.py $1 $GT ${CLASSES:+--classes $CLASSES} ${CODES:+--codes $CODES} --tag "$2" 2>&1 | grep -v Warning; }
sc seed_states.txt "board seeds"; sc refined_states.txt "discovery leaves"
# ---- (2) leaf tree, (3) merge-up
[ -s smerge/merge_tree.txt ] || $PY $W/state_merge.py smerge --rdir $RD --labels refined_states.txt --gt $GT --unfil --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "^== \|leaf pair" | cut -c1-200
log "merge-up along the leaf tree"; TAG=$TAG LEAVES=$O/refined_states.txt TREE=$O/smerge/merge_tree.txt PROJ=$RD/p.simple OUT=$O/bu GT=$GT PAR=${PAR:-8} bash $W/bu_merge.sbatch 2>&1 | grep -v "^\[.*eo:" | tail -40
log DONE
