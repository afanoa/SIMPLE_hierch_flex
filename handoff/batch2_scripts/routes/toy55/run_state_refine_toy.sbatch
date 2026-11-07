#!/bin/bash
#SBATCH -J toy10_state_refine
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The handoff test: an UNSUPERVISED K-state partition of the klem latent (weighted k-means on the posterior means) seeds
# SIMPLE's multi-state refinement (refine3D_states, local poses, 20->8 A); the refined state labels are scored against the
# truth per conformer (purity / recall, Hungarian). env: KDIR (latent dir), K (10), MAXITS (8), NP (10), RTAG (output dir
# suffix), STATEMODE (cells = the argmax cells as states), PROJSRC (project to start from, default the consensus refine3D)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464 OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; T=$W/sim; K=${K:-10}; NP=${NP:-10}; KDIR=${KDIR:-$W/lat/klem_junk}
D=$W/lat/state_refine_k$K${STATEMODE:+_$STATEMODE}${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') sr] $*"; }
RP=$(ls -td $W/proj/*_refine3D | head -1); PROJSRC=${PROJSRC:-$(ls $RP/*.simple | head -1)}; NPTCL=$(grep -v '^#' $T/gt_particles.txt | grep -c .)
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$KDIR/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z = L[:, 1:3]; gm = L[:, 5]; am = L[:, 3].astype(int); ok = gm > 0.5; rng = np.random.default_rng(0); X = z[ok]; K = $K
if '${STATEMODE:-}' == 'cells':   # the model's own states: argmax cell of every on-grid particle
    okc = ok & (am <= K); states = np.zeros($NPTCL, int); states[pind[okc] - 1] = am[okc]; np.savetxt('$D/states.txt', states, fmt='%d'); print('states = argmax cells: %d particles over %d states' % (okc.sum(), len(np.unique(am[okc])))); raise SystemExit
best = None
for seed in range(10):
    C = X[rng.choice(len(X), K, replace=False)]
    for it in range(100):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2); a = d.argmin(1); Cn = np.array([X[a == k].mean(0) if (a == k).any() else C[k] for k in range(K)])
        if np.allclose(Cn, C): break
        C = Cn
    cost = d[np.arange(len(X)), a].sum()
    if best is None or cost < best[0]: best = (cost, a.copy())
states = np.zeros($NPTCL, int); states[pind[ok] - 1] = best[1] + 1; np.savetxt('$D/states.txt', states, fmt='%d')
print('k-means K=%d on the posterior means: sizes %s; %d particles unassigned (free components)' % (K, np.bincount(best[1]), (~ok).sum()))
PYEOF
log "states file rows: $(wc -l < $D/states.txt) (project particles $NPTCL)"; cp -f $PROJSRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; tail -3 sel.log; exit 1; }
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/seed_ptcl3D.txt
$PY $W/state_purity.py $D/seed_ptcl3D.txt $T/gt_particles.txt --tag "SEED (klem k-means, before refinement)"
log "refine3D_states nstates=$K local 20->8 A maxits ${MAXITS:-8}"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-8} mskdiam=200 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt
$PY $W/state_purity.py $D/refined_ptcl3D.txt $T/gt_particles.txt --tag "REFINED (refine3D_states, hard labels)"
log "poses: $($PY $W/pose_error.py $D/refined_ptcl3D.txt $T/all.txt --tag poses)"
log DONE
