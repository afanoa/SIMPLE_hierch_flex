#!/bin/bash
#SBATCH -J oracle_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# oracle_refine = node_refine with the state labels given (env STATES: one label per row, 0 = left out, e.g. the true
# conformers) instead of k-means on the node board; K = number of labels. node_refine: per-node pose refinement in ONE half with the f692e family recipe (f692_family.sh steps 2-3): k-means K on the node
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
[ -s ${STATES:?} ] || $PY /mnt/beegfs/elmlund/afan/routes/toy55/oracle_states.py ${GT:?} ${SUBSET:?} $STATES || exit 1
n=$(awk '$1>0' $STATES | grep -c .); K=$(awk '$1>0' $STATES | sort -u | grep -c .); np_=$(( n / 8000 )); [ $np_ -lt 2 ] && np_=2; [ $np_ -gt 10 ] && np_=10; NP=${NP:-$np_}
NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/meta_toy.npz')['flags']))"); D=$OUT/refine; mkdir -p $D; cd $D; log(){ echo "[$(date '+%F %T') node_refine $ND] $*"; }
fallback(){ log "FAILED ($1): the node keeps its input project and poses"; ln -sfn $PROJ $OUT/proj.simple; ln -sfn $POSES $OUT/poses_full.txt; touch $D/FAILED; exit 0; }
[ -f $D/DONE ] && exit 0; mkdir -p $ND
cp $STATES $D/states.txt || fallback states; log "states from $STATES: $(awk '$1>0' $D/states.txt | sort -n | uniq -c | tr '\n' ' ')"
cp -f $PROJ p.simple || fallback copy
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || fallback selection
export OMP_NUM_THREADS=8; log "$n particles: refine3D_states nstates=$K local $LPS->8 A maxits ${MAXITS:-6} nparts $NP x 8 threads (${REFINE_EXTRA:-ml_reg=no})"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=$LPS lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 ${REFINE_EXTRA:-ml_reg=no} > refine.log 2>&1 && grep -a -q ">>> ITERATION" refine.log || fallback "refine3D_states: $(grep -a -m1 'ERROR\|STOP\|SIGSEGV' refine.log | cut -c1-120)"
RD=$(ls -td $D/*_refine3D_states | head -1); [ -s $RD/p.simple ] || fallback "no refined project"
rm -f $OUT/poses_full.txt   # a failed earlier attempt left a symlink to the INPUT table: never write through it
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
