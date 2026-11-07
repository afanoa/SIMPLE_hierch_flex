#!/bin/bash
#SBATCH -J node_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Per-node pose refinement in ONE half with the f692e family recipe (f692_family.sh steps 2-3): k-means K on the node
# board's posterior means (on-grid particles) -> refine3D_states K states, pose_policy=local, lpstart LPS -> lpstop 8,
# maxits 6, mskdiam 220, REFINE_EXTRA (ml_reg=no as in f692e). Only this half's particles and maps: the halves stay
# independent. Writes OUT/proj.simple and OUT/poses_full.txt (every row of the half); on failure the node keeps its
# input project and poses. env R (half run dir), PROJ, POSES (input project / pose table), ND (node dir: subset.txt,
# board/klem_latent.txt), OUT ($ND/refined), LPS (20), REFK (states, 4-10 by size), NP (2-10 by size), MAXITS (6), REFINE_EXTRA (ml_reg=no)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; R=${R:?}; PROJ=${PROJ:?}; POSES=${POSES:?}; ND=${ND:?}; OUT=${OUT:-$ND/refined}; LPS=${LPS:-20}
n=$(grep -c . $ND/subset.txt); k_=$(( n / 4000 )); [ $k_ -lt 4 ] && k_=4; [ $k_ -gt 10 ] && k_=10; K=${REFK:-$k_}; np_=$(( n / 8000 )); [ $np_ -lt 2 ] && np_=2; [ $np_ -gt 10 ] && np_=10; NP=${NP:-$np_}
NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/meta_toy.npz')['flags']))"); D=$OUT/refine; mkdir -p $D; cd $D; log(){ echo "[$(date '+%F %T') node_refine $ND] $*"; }
fallback(){ log "FAILED ($1): the node keeps its input project and poses"; ln -sfn $PROJ $OUT/proj.simple; ln -sfn $POSES $OUT/poses_full.txt; touch $D/FAILED; exit 0; }
[ -f $D/DONE ] && exit 0
$PY - <<PYEOF || fallback k-means
import numpy as np
L = np.loadtxt('$ND/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z = L[:, 1:3]; gm = L[:, 5]; ok = gm > 0.5; rng = np.random.default_rng(0); X = z[ok]; K = $K; best = None
for seed in range(10):
    C = X[rng.choice(len(X), K, replace=False)]
    for it in range(100):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2); a = d.argmin(1); Cn = np.array([X[a == k].mean(0) if (a == k).any() else C[k] for k in range(K)])
        if np.allclose(Cn, C): break
        C = Cn
    cost = d[np.arange(len(X)), a].sum()
    if best is None or cost < best[0]: best = (cost, a.copy())
st = np.zeros($NROWS, int); st[pind[ok] - 1] = best[1] + 1; np.savetxt('$D/states.txt', st, fmt='%d'); print('k-means K=%d on the node board posterior means (%d on-grid of %d): sizes %s' % (K, ok.sum(), len(ok), np.bincount(best[1])))
PYEOF
cp -f $PROJ p.simple || fallback copy
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || fallback selection
export OMP_NUM_THREADS=8; log "$n particles: refine3D_states nstates=$K local $LPS->8 A maxits ${MAXITS:-6} nparts $NP x 8 threads (${REFINE_EXTRA:-ml_reg=no})"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=$LPS lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 ${REFINE_EXTRA:-ml_reg=no} > refine.log 2>&1 && grep -a -q ">>> ITERATION" refine.log || fallback "refine3D_states: $(grep -a -m1 'ERROR\|STOP\|SIGSEGV' refine.log | cut -c1-120)"
RD=$(ls -td $D/*_refine3D_states | head -1); [ -s $RD/p.simple ] || fallback "no refined project"
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $OUT/poses_full.txt; [ $(grep -c . $OUT/poses_full.txt) -ge $NROWS ] || fallback "pose table"
ln -sfn $RD/p.simple $OUT/proj.simple; rm -f $RD/*part*.mrc 2>/dev/null
$PY - <<PYEOF 2>/dev/null
import numpy as np
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
st = np.loadtxt('$D/states.txt', dtype=int); want = st > 0
def load(f): return np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'x', 'y')] for i, l in enumerate(open(f)) if i < len(want) and want[i]])
A, Bv = load('$POSES'), load('$OUT/poses_full.txt')
def v(e): e1, e2 = np.radians(e[:, 0]), np.radians(e[:, 1]); return np.c_[np.cos(e1) * np.sin(e2), np.sin(e1) * np.sin(e2), np.cos(e2)]
ang = np.degrees(np.arccos(np.clip((v(A) * v(Bv)).sum(1), -1, 1))); sh = np.hypot(A[:, 2] - Bv[:, 2], A[:, 3] - Bv[:, 3])
print('[refined] %d particles: viewing-direction change median %.1f, p90 %.1f deg; shift change median %.2f px' % (len(A), np.median(ang), np.percentile(ang, 90), np.median(sh)))
PYEOF
touch $D/DONE; log DONE
