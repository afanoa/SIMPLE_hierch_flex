#!/bin/bash
#SBATCH -J toy10_r3s
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Standalone refine3D_states (K states from the coarse-window labels already in $OUT/r3s/p.simple), truth-scored. env: TAG, SIMDIR, CSV, K, NP
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/toyexp/pv; TAG=${TAG:-}; K=${K:-4}; NP=${NP:-10}; SIMDIR=${SIMDIR:-sim}
OUT=$W/lat$TAG; T=$W/$SIMDIR; GT=$T/gt_particles.txt; CSV=$W/${CSV:-rnasep_conformers10.csv}; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $CSV)
TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc $W/maps160/variants/conf$(printf %02d $c)_*.mrc 2>/dev/null | head -1; done | tr "\n" , | sed 's/,$//')
log(){ echo "[$(date '+%F %T') r3s$TAG] $*"; }
score(){ $PY - $1 $GT <<'PYEOF'
import sys, numpy as np
from itertools import permutations
lab = {}
for i, l in enumerate(open(sys.argv[1])):
    kv = dict(t.split('=') for t in l.split() if '=' in t); s = int(float(kv.get('state', 0)))
    if s > 0: lab[i + 1] = s
gt = {}; closed = {}
for l in open(sys.argv[2]):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1]); closed[int(w[0])] = int(w[2])
p = np.array(sorted(lab)); L = np.array([lab[i] for i in p]); t = np.array([gt[i] for i in p]); c = np.array([closed[i] for i in p]); ks = np.unique(L); keep = t > 0
ts = np.unique(t[keep]); Lk, tk = L[keep], t[keep]
best = max(((sum(np.sum((Lk == k) & (tk == s)) for k, s in zip(pm, ts[:len(pm)])), pm) for pm in permutations(ks, min(len(ks), len(ts)))), key=lambda x: x[0]) if len(ks) <= len(ts) else (0, ())
pur = sum(np.bincount(tk[Lk == k]).max() for k in ks if np.any(Lk == k)) / max(keep.sum(), 1)
cl = [(np.sum((L == k) & (c == 1)), k) for k in ks]; nc, kc = max(cl)
print('labelled %d (junk %.3f) | Hungarian accuracy vs %d conformers (K=%d) %.3f | weighted purity %.3f | closed: label %d holds %.2f at purity %.2f | sizes %s' %
      (len(p), (t == 0).mean(), len(ts), len(ks), best[0] / max(keep.sum(), 1), pur, kc, nc / max((c == 1).sum(), 1), nc / max(np.sum(L == kc), 1), ' '.join('%d:%d' % (k, np.sum(L == k)) for k in ks)))
PYEOF
}
D=$OUT/r3s; cd $D || exit 1; log "START $D"
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > seed_ptcl3D.txt; log "seed: $(score seed_ptcl3D.txt)"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=200 pgrp=c1 nparts=$NP nthr=8 > refine2.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR' refine2.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/ptcl3D_final.txt
log "refine3D_states (K=$K, local poses, 20->8 A): $(score $D/ptcl3D_final.txt)"
MS=$(ls $RD/recvol_state0?.mrc 2>/dev/null | tr "\n" , | sed 's/,$//'); [ -n "$MS" ] && $PY $P/pv_readout.py maps 1.464 15 200 $TRUTH $MS 2>&1 | grep -v " - " | head -6 | cut -c1-300
$PY $W/pose_error.py $D/ptcl3D_final.txt $T/all.txt --tag "poses after refine3D_states"
log DONE
